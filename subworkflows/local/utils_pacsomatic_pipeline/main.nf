/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    UTILITY FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Check non-empty parameters
def checkParameters() {
    def checkParamList = [
        [value: params.fasta, name: 'fasta']
    ]

    checkParamList.each { param ->
        if (!param.value) {
            log.error "Missing required parameter: --${param.name}"
            exit 1
        }
    }
}


// Check path parameters
def checkPathParameters() {
    def checkPathParamList = [params.fasta]
    checkPathParamList.each { param ->
        if (param) { file(param, checkIfExists: true) }
    }
}
