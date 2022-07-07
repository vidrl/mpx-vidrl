// Helper function to check if file exists and 
// return the Nextflow file(...) for staging
def check_file(file_path) {
    def fp = new File(file_path)
    assert fp.exists() : "File not found: ${fp}"
    return file(file_path)
}