scripts=".."
root=`dirname $scripts`
mkdir $root/software/build
automake_dir=/usr/local/Cellar/automake/1.15/share/automake-1.15

# Install modified ART
mkdir ../../software/build/art &&
    cd ../../software/build/art &&
    cp ../../art_illumina_src151.tar.gz . &&
    tar -xzf art_illumina_src151.tar.gz &&
    cp ../../art_illumina_src151-adapter-enabled.tar.gz . &&
    tar -xzf art_illumina_src151-adapter-enabled.tar.gz &&
    cd art_illumina_dir &&
    for f in config.sub config.guess install-sh depcomp missing INSTALL
    do
      rm $f
      ln -s $automake_dir/$f .
    done &&
    ./configure --prefix $root/software &&
    make &&
    make install &&
    cd ../../../scripts/simulation
