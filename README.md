## Title
Bioinformatic Toolkit
 
## Team Members
Alex Olson
Sola Yun
Will Hutcheon
 
## Nickname
Biokit

## Description (about 100 words)
This app will work as a tool kit to run common bioinformatic algorithms, using nucleotide string data as input. The first algorithm for the prototype will be creating and displaying a phylogenetic tree with persistent data storage.
As of now, the plan is to use the https://mygene.info/v3 api, to get the biological data.
This app will have 3 user classes, an employee/researcher, manager, and admin.
The employee will be able to upload data, and generate a report, as well as print out a report that they have previously generated (but not another employees), the manager will be able to print out any employees reports, and also tweak the customization settings for the report generation itself, along with add/delete/modify privileges. The admin will have control over adding/removing new users, changing passwords, and setting permissions.
The program will have a login/authentication system, along with persistent data storage.



> What are the goals of the app?
To have a collection of easy to use genetic tools for researchers to compile data and print reports.

> What problem will the app solve?


Handling nucleotide data can be very unwieldy, especially for newcomers, this app will abstract some of that complexity and generate data reports along with professional looking graphs quickly and efficiently. 

> What form will the app take to fulfill its goals?
As of now, it will be a web app created through Python’s Flask framework, along with HTML. We will initially have persistent data storage using a text file and CSV's, but will expand on that in later iterations 

> What sorts of features will the app have?
The app will be able to process multiple nucleotide strings and output a report with various information, including GC content, isochore data, the compilated phylogenetic tree, a dot matrix chart between two nucleotide strings, and others.
> 
![Use_case](https://github.com/aolson078/biokit/assets/69769089/97b19bdb-c369-4d3a-a883-24ca40f4b959)
