nextflow.enable.dsl = 2

params.carveme_container = "${projectDir}/../containers/carveme-1.6.6.sif"

process TEST_CARVEME_CONTAINER {

    container params.carveme_container

    output:
    path 'container-test.txt'

    script:
    """
    set -e

    carve --help > carve-help.txt

    python -c "import carveme"
    python -c "import reframed"
    python -c "import pyscipopt"
    python -c "import libsbml"

    {
        echo "Container: ${params.carveme_container}"
        python -c "from importlib.metadata import version; print('CarveMe:', version('carveme'))"
        python -c "from importlib.metadata import version; print('ReFramed:', version('reframed'))"
        diamond version
    } > container-test.txt
    """
}

workflow {
    TEST_CARVEME_CONTAINER()
}
