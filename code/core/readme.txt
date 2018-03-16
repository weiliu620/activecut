To compile the code:
1) run: export BOOST_ROOT=/home/sci/weiliu/packages/boost_1_53_0/build/wukong
We need this because a bug of cmake 2.8.8 can not find_package(boost)

2) For in-source build, in current directory, run: ccmake ..
   Give the itk binary build directory in ccmake window. If you don't want to build itk, you can use the one in my directory: /home/sci/weiliu/packages/itk_4.3.1_bin/wukong 

 The boost dir is still 'not found', but it doesn't matter since we did step 1)

3) 'c' to configue and 'g' to generate makefile.

4) make

============================================ 

grabcutseg take a 4d iamge as input. The various channels of images can be
easily merged into one single 4D image by using some other tools (such as
convertITKformat). It is conceptually easy for grabcutseg to use 4D image
input. Ineed the algorithm does not care if the i'th image is T1 or T2. All
images are treated equally for now. The output will be 4D images, too.



All index and labels inside the code is 0-based. All label maps saved in file
are 1-based. Therefore. The code takes care of the conversion between them. For
example, when convert the gmm_labels into a image file. The 0-based labels with
range in [0, K-1] is converted to 1-based labels with rage in [1, K].
