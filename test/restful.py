# testing for RESTful API
import os.path
import requests
import datetime
import pprint

path = "/media/pika/Workbench/workspace/wgMLST/01--ContigAnnotation/testInputFiles/Assembly"
timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')

files = {'file': ('Assembly_1.fa', open(os.path.join(path, 'Assembly_1.fa'), 'rb'), 'text/x-fasta')}
form = {"created": timestamp, "filename": "test1"}
r = requests.post("http://127.0.0.1:5000/api/uploads", data=form, files=files)
# r = requests.post("http://10.0.30.80:5000/api/uploads", data=form, files=files)
r
pprint.pprint(r.json())

files = {'file': ('Assembly_2.fa', open(os.path.join(path, 'Assembly_2.fa'), 'rb'), 'text/x-fasta')}
form = {"created": timestamp, "filename": "test2"}
r = requests.post("http://127.0.0.1:5000/api/uploads", data=form, files=files)
r
pprint.pprint(r.json())
