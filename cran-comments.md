Uwe Ligges <ligges@statistik.tu-dortmund.de>
Fri 1/10/2020 7:58 PM
You; CRAN?

Thanks, we see lots of problems.

First of all note, an example must not take longer than 5 sec.
The whole package check must not tak longer than 10 min.

Moire details follow:

Flavor: r-devel-linux-x86_64-debian-gcc, r-devel-windows-ix86+x86_64
Check: CRAN incoming feasibility, Result: NOTE
?? Maintainer: 'Weimin Zhang <zhangwm@hotmail.com>'

?? New submission

?? Version contains large components (0.0.0.9000)

?? Possibly mis-spelled words in DESCRIPTION:
???? VUS (3:170, 15:176)

?? Package has a VignetteBuilder field but no prebuilt vignette index.

?? The build time stamp is missing.

Flavor: r-devel-linux-x86_64-debian-gcc, r-devel-windows-ix86+x86_64
Check: for hidden files and directories, Result: NOTE
?? Found the following hidden files and directories:
???? .Rhistory
???? .gitignore
???? vignettes/.gitignore
???? .Rproj.user
???? .git
?? These were most likely included in error. See section 'Package
?? structure' in the 'Writing R Extensions' manual.

?? CRAN-pack does not know about
???? .git

Flavor: r-devel-linux-x86_64-debian-gcc, r-devel-windows-ix86+x86_64
Check: DESCRIPTION meta-information, Result: NOTE
?? Checking should be performed on sources prepared by 'R CMD build'.

Flavor: r-devel-linux-x86_64-debian-gcc, r-devel-windows-ix86+x86_64
Check: top-level files, Result: NOTE
?? Non-standard files/directories found at top level:
???? 'Read-and-delete-me' 'sorocs.Rproj'

Flavor: r-devel-linux-x86_64-debian-gcc, r-devel-windows-ix86+x86_64
Check: package subdirectories, Result: WARNING
?? Found the following directory with the name of a version control 
directory:
???? ./.git
?? These should not be in a package tarball.

Flavor: r-devel-linux-x86_64-debian-gcc, r-devel-windows-ix86+x86_64
Check: R code for possible problems, Result: NOTE
?? sorocs: no visible binding for '<<-' assignment to 'dnormfun'
?? sorocs: no visible binding for '<<-' assignment to 'pnormfun'

Flavor: r-devel-windows-ix86+x86_64
Check: files in 'vignettes', Result: WARNING
?? Files in the 'vignettes' directory but no files in 'inst/doc':
???? 'sorocs-vignette.Rmd'

Flavor: r-devel-windows-ix86+x86_64
Check: running examples for arch 'i386', Result: NOTE
?? Examples with CPU (user + system) or elapsed time > 10s
??????????? user system elapsed
?? sorocs 759.01?? 54.3? 817.27

Flavor: r-devel-windows-ix86+x86_64
Check: running examples for arch 'x64', Result: FAIL
?? Check process probably crashed or hung up for 20 minutes ... killed
?? Most likely this happened in the example checks (?),
?? if not, ignore the following last lines of example output:

?? ======== End of example output (where/before crash/hang up occured ?) 
========

Flavor: r-devel-linux-x86_64-debian-gcc
Check: files in ‘vignettes’, Result: WARNING
?? Files in the 'vignettes' directory but no files in 'inst/doc':
???? 'sorocs-vignette.Rmd'

Flavor: r-devel-linux-x86_64-debian-gcc
Check: examples, Result: NOTE
?? Examples with CPU (user + system) or elapsed time > 10s
??????????? user system elapsed
?? sorocs 262.76? 8.344 271.208

Flavor: r-devel-linux-x86_64-debian-gcc
Check: package vignettes, Result: WARNING
?? dir.exists(dir) is not TRUE
?? Package vignette without corresponding single PDF/HTML:
????? 'sorocs-vignette.Rmd'

Please fix and resubmit.

Best,
Uwe Ligges

___________________________________________________________________________
ligges@statistik.tu-dortmund.de
Wed 1/22/2020 2:01 PM
You; CRAN-submissions@R-project.org
"00details.log".log
847 bytes
Dear maintainer,
 Flavor: r-devel-linux-x86_64-debian-gcc, r-devel-windows-ix86+x86_64

Check: CRAN incoming feasibility, Result: NOTE

  Maintainer: 'Weimin Zhang <zhangwm@hotmail.com>'

  

  New submission

  

  Possibly mis-spelled words in DESCRIPTION:

    VUS (5:46, 17:176)



Flavor: r-devel-linux-x86_64-debian-gcc, r-devel-windows-ix86+x86_64

Check: R code for possible problems, Result: NOTE

  sorocs: no visible binding for '<<-' assignment to 'dnormfun'

  sorocs: no visible binding for '<<-' assignment to 'pnormfun'
  
#
#search "no visible binding for '<<-' assignment to" got:
#https://community.rstudio.com/t/how-to-solve-no-visible-binding-for-global-variable-note/28887
#you can use the .data pronoun cf this dplyr vignette , " If this function is in a package, using .data also prevents #R CMD check from giving a NOTE about undefined global variables (provided that you’ve also imported rlang::.data #with @importFrom rlang .data )."
#  rlang::.data 
  
https://cran.r-project.org/web/packages/future/vignettes/future-4-issues.html
https://cran.r-project.org/doc/manuals/R-lang.html#Scope-of-variables
http://adv-r.had.co.nz/Functions.html
https://stackoverflow.com/questions/12227689/function-within-function-in-r
From: Max Turgeon <Max.Turgeon@umanitoba.ca>
Sent: Thursday, February 20, 2020 10:50 AM
To: Weimin Zhang <zhangwm@hotmail.com>; r-package-devel@r-project.org <r-package-devel@r-project.org>
Subject: Re: no visible binding for '<<-' assignment note 
 
Hi Weimin,

From the R docs, if there is no visible binding for the deep assignment, meaning there is the variable "dnormfun" not already defined in the parent environment, then an assignment is made in the global assignment. Is this really what you want to do, define a function "dnormfun" that is available in the user's global environment? If yes, then you'll need to make a case for it in your CRAN submision comments. If no, then you can probably just change "<<-" to "<-".

Best,


Max Turgeon
Assistant Professor
Department of Statistics
Department of Computer Science
University of Manitoba
maxturgeon.ca
++++++++++++++++++++++++++++++++++++++++++++++++++++++++





Max Turgeon <Max.Turgeon@umanitoba.ca>
Fri 2/21/2020 1:56 PM
You; r-package-devel@r-project.org

The only thing I know about your package is the code you've shared so far, so I may be wrong. If I understand correctly, you're using kronecker and dnormfun to evalute the normal density at several values of mu and sigma, but a fixed point x. And you're fixing that point x in dnormfun, with the possibility of redefining dnormfun if you need a new point x. If that is the case, I think you can simply do something like

kronecker(X, Y, FUN = dnorm, x = gridY[i]).

But again, maybe I'm missing something.


Max Turgeon

