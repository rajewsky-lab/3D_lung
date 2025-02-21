#Please adjust STIM and file paths 
mkdir ../stim_files/

## Resave input files in n5 format for efficient storage and access
st-resave \
-i ../stim_files/section_4_locations.csv,../stim_files/section_4_reads.csv,../stim_files/section_4_celltypes.csv,section_4.n5 \
-i ../stim_files/section_10_locations.csv,../stim_files/section_10_reads.csv,../stim_files/section_10_celltypes.csv,section_10.n5 \
-i ../stim_files/section_16_locations.csv,../stim_files/section_16_reads.csv,../stim_files/section_16_celltypes.csv,section_16.n5 \
-i ../stim_files/section_22_locations.csv,../stim_files/section_22_reads.csv,../stim_files/section_22_celltypes.csv,section_22.n5 \
-i ../stim_files/section_28_locations.csv,../stim_files/section_28_reads.csv,../stim_files/section_28_celltypes.csv,section_28.n5 \
-i ../stim_files/section_34_locations.csv,../stim_files/section_34_reads.csv,../stim_files/section_34_celltypes.csv,section_34.n5 \
-o cosmx.n5

## Align seach section with the 2 consecutive sections
st-align-pairs -i cosmx.n5 -n=100 -r 2  --overwrite -d 'section_4,section_16,section_22,section_28,section_34'

## Find optimal global alignment
st-align-global -i cosmx.n5 --absoluteThreshold 100 -sf 0.5 --lambda 0.5 --skipICP  -d 'section_4,section_16,section_22,section_28,section_34'