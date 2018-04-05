# make sure to push dev first
git checkout dev
git add --all .
git commit -m "small change"
git push origin dev

# cleanup current tmp
rm -fR /tmp/primerDAFT-tmp
mkdir /tmp/primerDAFT-tmp

# copy all nessesary files to tmp
cp {.gitignore,LICENSE,MANIFEST.in,README.md,make_docs.sh,requirements.txt,setup.py} /tmp/primerDAFT-tmp/
cp -R bin /tmp/primerDAFT-tmp/
cp -R docs /tmp/primerDAFT-tmp/
cp -R primerDAFT /tmp/primerDAFT-tmp/
cp -R sphinx /tmp/primerDAFT-tmp/
cp -R test /tmp/primerDAFT-tmp/

# go to master branch
git checkout master
git pull

# remove nessasary files
rm .gitignore LICENSE MANIFEST.in README.md make_docs.sh requirements.txt setup.py
rm -R bin
rm -R docs
rm -R primerDAFT
rm -R sphinx
rm -R test

# copy back all files in tmp to master and push
cp -a /tmp/primerDAFT-tmp/. ./
git add --all .
git commit -m "$1"
git push origin master

# come back to dev
git checkout dev
