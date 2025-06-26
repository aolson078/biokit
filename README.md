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
As of now, the plan is to use BioPython's Entrez library to interact with the NCBI (National Center for Biotechnology Information) API, to get the biological data.
This app will have 3 user classes, an employee/researcher, manager, and admin.
The employee will be able to upload data, and generate a report, as well as print out a report that they have previously generated (but not other employees), the manager will be able to print out any employee's reports, and also tweak the customization settings for the report generation itself, along with add/delete/modify privileges. The admin will have control over adding/removing new users, changing passwords, and setting permissions.
The program will have a login/authentication system and persistent data storage.



> What are the goals of the app?
To have a collection of easy-to-use genetic tools for researchers to compile data and print reports.

> What problem will the app solve?


Handling nucleotide data can be very unwieldy, especially for newcomers, this app will abstract some of that complexity and generate data reports along with professional-looking graphs quickly and efficiently. 

> What form will the app take to fulfill its goals?
Currently, it will be a web app created through Python’s Flask framework, along with HTML. We will initially have persistent data storage using a text file and CSV's, but will expand on that in later iterations 

> What sorts of features will the app have?
The app will be able to process multiple nucleotide strings and output a report with various information, including GC content, isochore data, the compilated phylogenetic tree, a dot matrix chart between two nucleotide strings, and others.
> 
![Use_case](https://github.com/aolson078/biokit/assets/69769089/97b19bdb-c369-4d3a-a883-24ca40f4b959)



Default login:
email: 1@gmail.com
password: 11111

To run:
Run file "routes.py"
Click ip address in python console (127.0.0.1:5000)
Login in with own or default user.

## Background Tasks
Heavy computations like report compilation run as asynchronous Celery jobs. Start a worker with `celery -A tasks.celery worker --loglevel=info` and use `/compile_report_async` to queue jobs. Check progress via `/task_status/<task_id>`.

## API Usage
Programmatic access is available through JSON endpoints:
- `POST /api/create_record` to add a record
- `POST /compile_report_async` to generate a report asynchronously
- `GET /api/report/<id>` to retrieve report data

## Docker
Build the application with Docker and run using docker-compose:
```bash
docker-compose up --build
```
This starts the web service, Celery worker, and Redis broker.
