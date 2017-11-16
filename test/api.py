import requests
import datetime
import os
from profiler.models import File

batch_id = '6af969dc5a604ca695410d7de8ccc1c6'
batch_created = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
ip = "localhost"
port = 5000


upload_files = list(File.objects.filter(batch_id__exact=batch_id))
url = "http://{0}:{1}/api/uploads".format(ip, port)
for f in upload_files:
    form = {'batch_id': batch_id, 'filename': os.path.splitext(os.path.basename(f.file.path)), 'created': batch_created}
    files = [('file', open(f.file.path, 'rb'))]
    requests.post(url, data=form, files=files)
url = "http://127.0.0.1:5000/api/profiling/{}".format(batch_id)
requests.get(url)
