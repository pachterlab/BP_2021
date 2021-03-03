wget -O ../../reference.tar.gz --content-disposition https://data.caltech.edu/tindfiles/serve/846e5231-109e-4672-a785-9496726c9ecf

tar -xvf ../../reference.tar.gz -C ../../

rm ../../reference.tar.gz

gunzip -r ../../reference/
