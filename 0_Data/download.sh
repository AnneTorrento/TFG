#!/bin/bash

confirm() {
    # call with a prompt string or use a default
    read -r -p "${1:-Are you sure? [y/N]} " response
    case "$response" in
        [yY][eE][sS]|[yY]) 
            true
            ;;
        *)
            false
            ;;
    esac
}

confirm "Download Atlas data? [y/N]" && curl https://content.cruk.cam.ac.uk/jmlab/atlas_data.tar.gz > atlas_data.tar.gz

# arrayexpress link:
# https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6967/E-MTAB-6967.processed.1.zip
