# Development Environment Setup

## Cloning the repository to your local computer
```shell
git clone git@github.com:briancohn/retina.git && cd retina && r
```

If you want to load a specific branch, run `git checkout -b iss13 && git pull origin iss13`

```r
install.packages('devtools'); library(devtools); install(); test(); load_all();
```

## How to commit a new change from the terminal if you are a collaborator on the project.
```shell
# Get an issue number (X) by github.com/bc/retina/issues/new
git checkout -b issX
git add R/some_file_you_edited.R
git commit -m '#13 I made a change to the file so it now does X'
git push origin issX
```
Be specific on the changes you made, & be sure to note what issue you are referring to via `#`.

# Want to help contribute?
Let me know at `brian.cohn@usc.edu` and I can get you up to speed so you can implement features and extend the codebase!
