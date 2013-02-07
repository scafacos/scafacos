# Assume we want to install them below $HOME/local.
myprefix=$HOME/local

# Show where the auto-tools will be installed to
echo "auto-tools will be installed to: $myprefix"

# Create a temporary directory to store the tar-balls
echo "creating temporary directory"
mkdir tmp
cd tmp

# Do the following in a scratch directory.
echo "downloading auto-tools"
wget http://ftp.gnu.org/gnu/m4/m4-1.4.16.tar.gz
wget http://ftp.gnu.org/gnu/autoconf/autoconf-2.68.tar.gz
wget http://ftp.gnu.org/gnu/automake/automake-1.11.3.tar.gz
wget http://ftp.gnu.org/gnu/libtool/libtool-2.4.2.tar.gz
echo "extracting auto-tools"
gzip -dc m4-1.4.16.tar.gz | tar xvf -
gzip -dc autoconf-2.68.tar.gz | tar xvf -
gzip -dc automake-1.11.3.tar.gz | tar xvf -
gzip -dc libtool-2.4.2.tar.gz | tar xvf -
echo "installing auto-tools to $myprefix"
cd m4-1.4.16
./configure -C --prefix=$myprefix && make && make install
export PATH=${myprefix}:${PATH}
cd ../autoconf-2.68
./configure -C --prefix=$myprefix && make && make install
cd ../automake-1.11.3
./configure -C --prefix=$myprefix && make && make install
cd ../libtool-2.4.2
./configure -C --prefix=$myprefix && make && make install

# Clean up the temporary directory and downloaded tar-balls
cd ../..
echo "deleting temporary directory `pwd`/tmp"
rm -rf tmp

# Notification to update the PATH variable
echo "$myprefix/bin was added to your PATH variable to use the newly installed auto-tools"
