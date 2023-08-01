#!/bin/bash
directories={
        ~/../../media/studentsgh129/Seagate\ Expansion\ Drive/ABBY.windows_surface_pro/FASTQ.files/AD/
        ~/../../media/studentsgh129/Seagate\ Expansion\ Drive/ABBY.windows_surface_pro/FASTQ.files/Old/
        ~/../../media/studentsgh129/Seagate\ Expansion\ Drive/ABBY.windows_surface_pro/FASTQ.files/Young/
}

for dir in "$directories"
do
        cd "$dir"

        for file in ./*.fastq.gz
        do
                salmon quant -i ~/salmon_V43_index/salmon_v43_index/ -l A -r $file --validateMappings --gcBias --seqBias
	done
done
