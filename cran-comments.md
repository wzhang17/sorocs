Uwe Ligges <ligges@statistik.tu-dortmund.de>
Fri 1/10/2020 7:58 PM
You; CRAN?

Thanks, we see lots of problems.

First of all note, an example must not take longer than 5 sec.
The whole package check must not tak longer than 10 min.

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

  Check: R code for possible problems, Result: NOTE

  sorocs: no visible binding for '<<-' assignment to 'dnormfun'

  sorocs: no visible binding for '<<-' assignment to 'pnormfun'
  
#
__________________________________________________________________

From: Max Turgeon <Max.Turgeon@umanitoba.ca>
Sent: Thursday, February 20, 2020 10:50 AM
To: Weimin Zhang <zhangwm@hotmail.com>; r-package-devel@r-project.org <r-package-devel@r-project.org>
Subject: Re: no visible binding for '<<-' assignment note 
 
Hi Weimin,

From the R docs, if there is no visible binding for the deep assignment, meaning there is the variable "dnormfun" not already defined in the parent environment, then an assignment is made in the global assignment. Is this really what you want to do, define a function "dnormfun" that is available in the user's global environment? If yes, then you'll need to make a case for it in your CRAN submision comments. If no, then you can probably just change "<<-" to "<-".


The only thing I know about your package is the code you've shared so far, so I may be wrong. If I understand correctly, you're using kronecker and dnormfun to evalute the normal density at several values of mu and sigma, but a fixed point x. And you're fixing that point x in dnormfun, with the possibility of redefining dnormfun if you need a new point x. If that is the case, I think you can simply do something like

kronecker(X, Y, FUN = dnorm, x = gridY[i]).

But again, maybe I'm missing something.


Max Turgeon
Assistant Professor
Department of Statistics
Department of Computer Science
University of Manitoba
maxturgeon.ca
++++++++++++++++++++++++++++++++++++++++++++++++++++++++

——————————————————————————————————————————————————————————————————————
Uwe Ligges <ligges@statistik.tu-dortmund.de>
Wed 2/26/2020 7:49 AM
You; CRAN

Thanks, we see:


   Found the following (possibly) invalid URLs:
     URL: doi:10.1111/biom.12997
       From: inst/doc/sorocs-vignette.html
       Message: Invalid URI scheme (use \doi for DOIs in Rd markup)

Within URL markup, you have to specify the full URL, of course.

Please fix and resubmit.

Best,
Uwe Ligges

__________________________________________________________________

Jelena Saf <jelena.saf@jku.at>
Wed 3/4/2020 2:22 PM
You; CRAN

Thanks,

Please limit your submission comments to a short summary of what you 
have changed since your last submission or why you decided not to change 
certain things that we recommended we do. We can't read all of this.

Please shorten the title to a maximum of 65 characters.
Acronyms can be used on their own in the title as long as they are 
explained in the description field.

Please do not start your description field in the DESCRIPTION file with 
phrases like 'This is a R package', 'This package', the package name or 
similar.

Only capitalize names and sentence beginnings in your description.
Also please have a space after periods: --> constraints. The

Please fix and resubmit, and document what was changed in the submission 
comments.

Best,
Jelena Saf