
// Helper function to check if file exists and 
// return the Nextflow file(...) for staging
def check_file(file_path) {
    def fp = new File(file_path)
    assert fp.exists() : "File not found: ${fp}"
    return file(file_path)
}

// Helper function to check if the primer scheme directory
// conforms to required formatting - this is the same check
// that is conducted in `minion.py` but fails before the 
// pipeline is run and does not trigger an obscure download
// without informative error message
def validate_primer_scheme(file_path) {

    def fp = new File(file_path)

    if (!fp.exists()) {
        println("Primer scheme directory not found: ${fp}")
        System.exit(1)
    }

    primer_scheme_dir = new File(fp.getAbsolutePath());
    primer_scheme_name_dir = primer_scheme_dir.getParentFile();
    primer_scheme_base_dir = primer_scheme_name_dir.getParentFile().getAbsolutePath();

    // println("Primer scheme directory: $primer_scheme_dir")
    // println("Parent directory: $primer_scheme_name_dir")
    // println("Primer scheme base directory: $primer_scheme_base_dir")

    primer_version = primer_scheme_dir.getName();
    scheme_name = primer_scheme_name_dir.getName();

    // println("Primer scheme version name: $primer_version")
    // println("Primer scheme name: $scheme_name")

    if (primer_version.startsWith("V")) {
        version = primer_version.minus("V")
    } else {
        println("Could not parse primer version: ${primer_version}")
        System.exit(1)
    }

    def fasta = new File(file_path, "${scheme_name}.reference.fasta")
    def bed = new File(file_path, "${scheme_name}.scheme.bed")

    if (!fasta.exists()) {
        println("Scheme reference file not found: ${fasta}")
        System.exit(1)
    }

    if (!bed.exists()) {
        println("Primer scheme file not found: ${bed}")
        System.exit(1)
    }

    // println("Primer scheme name: $scheme_name")
    // println("Primer scheme version: $version")

    return [tuple(file(primer_scheme_base_dir), "${scheme_name}/${primer_version}"), file(bed)]

}

def barcode_string_to_integer(range_string){
    try {
        return Integer.parseInt(range_string.replaceAll("\\s",""));  // strip whitespaces
    } catch (NumberFormatException e) {
        println("Could not convert barcode (value: $range_string) to integer.")
        System.exit(1)
    }
}

// Helper function to get a list of barcodes from user input
def get_barcodes_tails(barcode_str) {

    if (barcode_str instanceof Integer){
        println("Barcodes can only be specified as range or comma-delimited list.")
        System.exit(1)
    }

    if (barcode_str.contains("-")){
        // Barcodes are listed as range
        (start, stop) = barcode_str.split('-');

        barcode_start = barcode_string_to_integer(start);
        barcode_stop = barcode_string_to_integer(stop);

        if (barcode_stop <= barcode_start){
            println("End of range ($barcode_stop) cannot be equal or smaller than start of range ($barcode_start).")
            System.exit(1)
        }

        return (barcode_start..barcode_stop).collect { String.format('%02d', it) }

    } else if (barcode_str.contains(",")) {
        // Barcodes are listed as comma-delimited list
        return barcode_str.split(',').collect { 
            barcode = barcode_string_to_integer(it);
            String.format('%02d', barcode)
        }
    } else {
        println("Could not detect '-' for inclusive range or ',' for list of barcodes.")
        System.exit(1)
    }

}


// NB: empty warning is for some reason not 
// triggered when the channel is passed into
// processes after - this may be because the
// ifEmpty closure still emits / creates empty
// channel despite also exiting with error code

// For better readability in main workflow
def get_fastq_gather(file_glob){
    fastq_files = channel.fromPath(file_glob, type: 'file') | collect | map { tuple(params.fastq_id, it) }
    fastq_files | ifEmpty { exit 1, "Could not find any read files for processing" }
    return fastq_files
}

// For better readability in main workflow
def get_fastq_dirs(fastq_dir, fastq_ext, barcode_str){

    barcode_tails = null
    if (barcode_str){
        println("Filtering subdirectories by user-specified barcodes: ${barcode_str}")
        barcode_tails = get_barcodes_tails(barcode_str)
    }

    fastq_files = channel.fromPath("${fastq_dir}/*", type: 'dir') | map { dir -> 

        if (barcode_tails){
            found_barcode = barcode_tails.findAll { ext -> dir.baseName.endsWith(ext) }
            if (!found_barcode){
                return  // do not use this directory
            }
        }

        def files = []
        dir.eachFile { file -> 
            if (file.getName().endsWith("${fastq_ext}")){
                files << file
            }
        }
        if (files) {
            // println("Found sample directory: ${dir.baseName} (number of files: ${files.size()})")
            return tuple(dir.baseName, files) 
        }
    }

    fastq_files | ifEmpty { exit 1, "Could not find any suitable subdirectories containing read files for processing" }
    return fastq_files
}

// Sample sheet input

def get_samples(fastq_dir, fastq_ext, sample_sheet){

    fastq_files = channel.fromPath(sample_sheet) | splitCsv(header:true) | map { row -> 

        barcode_dir = file("${fastq_dir}/${row.barcode}")


        if (!barcode_dir.exists()){
            println("Barcode directory does not exist: ${barcode_dir} -> ${row.sample_id}")
            return 
        }

        def files = []
        barcode_dir.eachFile { file -> 
            if (file.getName().endsWith("${fastq_ext}")){
                files << file
            }
        }

        if (files) {
            // println("Found sample directory: ${row.barcode} -> ${row.sample_id} (number of files: ${files.size()})")
            return tuple(row.sample_id, files) 
        }
    
    }

    fastq_files | ifEmpty { exit 1, "Could not find any suitable subdirectories specified in sample sheet" }

    return fastq_files

}

// Wrapper for getting fastq files
def get_fastq_files(fastq_gather, fastq_id, fastq_dir, fastq_ext, barcodes, sample_sheet){


    if (fastq_gather) {
        println("Collecting read files: $fastq_gather")
        if (!fastq_id){
            println("Please provide an output identifier for the sample (--fastq_id)")
            System.exit(1)
        }
        fastq_files = get_fastq_gather(fastq_gather)
    } else if (fastq_dir) {
        println("Collecting read files with extension '$fastq_ext' in subdirectories of: $fastq_dir")
         if (!fastq_ext){
            println("Please provide a file extension to gather read files in subdirectories (--fastq_ext)")
            System.exit(1)
        }
        if (sample_sheet){
            println("Collecting read files from sample sheet: $sample_sheet")
            // Subset barcodes by sample sheet
            fastq_files = get_samples(fastq_dir, fastq_ext, sample_sheet)
        } else {
            // Collect all barcodes
            fastq_files = get_fastq_dirs(fastq_dir, fastq_ext, barcodes)
        }

    } else {
        println("Please provide read file inputs (--fastq_gather | --fastq_dir)")
        System.exit(1)
    }

}


// Sample sheet input

def get_samples_paired(fastq_dir, sample_sheet){

    fastq_files = channel.fromPath(sample_sheet) | splitCsv(header:true) | map { row -> 

        forward = file("${fastq_dir}/${row.forward}")
        reverse = file("${fastq_dir}/${row.reverse}")

        if (!sample_dir.exists()){
            println("Barcode directory does not exist: ${sample_dir} -> ${row.sample_id}")
            return 
        }

        if (!forward.exists()){
            println("Forward read file does not exist: ${forward}")
            return 
        }
        if (!reverse.exists()){
            println("Reverse read file does not exist: ${reverse}")
            return 
        }

        return tuple(row.sample_id, forward, reverse)
    
    }

    fastq_files | ifEmpty { exit 1, "Could not find any suitable subdirectories specified in sample sheet" }

    return fastq_files

}
