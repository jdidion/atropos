scripts=".."
root=`pwd`
mkdir $root/software/build
# assuming automake was installed using homebrew
automake_dir=/usr/local/Cellar/automake/1.15/share/automake-1.15

# Install modified ART
mkdir art/build &&
    cd art/build &&
    tar -xzf ../../art_illumina_src151.tar.gz -C . &&
    tar -xzf ../../art_illumina_src151-adapter-enabled.tar.gz -C . &&
    cd art_illumina_dir &&
    for f in config.sub config.guess install-sh depcomp missing INSTALL
    do
      rm $f
      ln -s $automake_dir/$f .
    done &&
    ./configure --prefix $root/art &&
    make &&
    make install
