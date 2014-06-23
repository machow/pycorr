# FIRST STEPS (can ignore if already have run script before)
# remove current gh-pages branch
git push origin :gh-pages
git checkout --orphan gh-pages
git rm *
git rm .*
echo "first commit" > index.html
git commit -am "first commit"
git push origin gh-pages

# UPDATING DOCS (run each time to push updates)
cd docs && rm -rf build/html
mkdir -p build/html
git clone --branch gh-pages git@github.com:machow/pycorr.git build/html
make html
cd build/html && touch .nojekyll
git commit -am "updated docs"
