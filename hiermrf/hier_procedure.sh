# Use SW to generate label image.
./swsim  --beta 0.9 -k 3 -m ./mask.nii.gz -g sw.nii.gz -b 5000 -n 20 -s -v 1 -a swsamples.nii.gz

# create noise scalar image.
~/projects/groupising/genobsimage -l sw.nii.gz -r 0.1 -o obs.nii.gz 

./singlemrf -b 20 -n 5 --emiter 30 --beta 0.9 -k 3 -s -i sw.nii.gz -f obs.nii.gz --samplefile swsamples.nii.gz -g swout.nii.gz -v 1


./singlemrf -b 20 -n 5 --emiter 30 --beta 0.9 -k 3  -i sw.nii.gz -f obs.nii.gz --samplefile gibbssamples.nii.gz -g gibbsout.nii.gz -v 1

~/projects/util/rain -r sw.nii.gz -m mask.nii.gz -i swout.nii.gz
~/projects/util/rain -r sw.nii.gz -m mask.nii.gz -i gibbsout.nii.gz

