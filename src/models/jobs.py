import os
import shutil
import sqlite3
import uuid
from datetime import datetime
from enum import Enum

from src.utils import files

TIME_FMT = "%Y-%m-%d %H:%M:%S"


def connect_db(db_dir):
    conn = sqlite3.connect(files.joinpath(db_dir, "data.db"))
    return conn


def create_table(db):
    db.execute("CREATE TABLE jobs (jobid text, create_time text, job_type text);")


class JobType(Enum):
    PGDB = {"str": "pgdb"}
    WGMLST = {"str": "wgmlst"}

    def to_str(self):
        return self.value["str"]

    @classmethod
    def all(cls):
        return [cls.PGDB, cls.WGMLST]

    @classmethod
    def parse(cls, s):
        for x in cls.all():
            if x.value["str"] == s:
                return x


class Job:
    def __init__(self, jobid, create_time, job_type):
        self._jobid = jobid
        self._create_time = create_time
        self._job_type = job_type

    @property
    def jobid(self)-> str:
        return self._jobid

    @property
    def create_time(self)-> str:
        return self._create_time.strftime(TIME_FMT)

    @property
    def job_type(self)-> str:
        return self._job_type.to_str()

    def to_tuple(self):
        return self.jobid, self.create_time, self.job_type


class JobManager:
    def __init__(self, db_dir):
        self.jobs = []
        self.db_dir = db_dir
        self._db = connect_db(db_dir)
        self._cur = self._db.cursor()
        self.load_jobs()

    def start_job(self, job_type: JobType):
        jobid = self._create_uuid()
        create_time = datetime.now()
        j = Job(jobid, create_time, job_type)
        self.jobs.append(j)
        job_dir = self._create_folder(j.job_type, j.jobid, self.db_dir)
        self._append2db(self._cur, j)
        self._db.commit()
        return jobid, job_dir

    def load_jobs(self):
        self._cur.execute("SELECT * FROM jobs;")
        job_lst = self._cur.fetchall()
        if job_lst:
            for j in job_lst:
                jobid, create_time, job_type = j
                create_time = datetime.strptime(create_time, TIME_FMT)
                job_type = JobType.parse(job_type)
                self.jobs.append(Job(jobid, create_time, job_type))

    def get_jobs(self, job_type: JobType):
        return [j for j in self.jobs if j.job_type == job_type.to_str()]

    def delete_job(self, job: Job=None, jobid: str=None):
        if job:
            pass
        elif jobid:
            for j in self.jobs:
                if j.jobid == jobid:
                    job = j
                    break
        else:
            raise Exception("Must assign one of job or job id for deletion.")

        self.jobs.remove(job)
        self._delete_folder(job.job_type, job.jobid, self.db_dir)
        self._delete_from_db(self._cur, job.jobid)
        self._db.commit()

    def close(self):
        self._db.commit()
        self._db.close()

    @classmethod
    def _create_folder(cls, job_type: str, jobid: str, db_dir):
        job_dir = files.joinpath(db_dir, job_type, jobid)
        os.mkdir(job_dir)
        return job_dir

    @classmethod
    def _delete_folder(cls, job_type: str, jobid: str, db_dir):
        shutil.rmtree(files.joinpath(db_dir, job_type, jobid))

    @classmethod
    def _append2db(cls, db, job: Job):
        db.execute("INSERT INTO jobs VALUES (?,?,?)", job.to_tuple())

    @classmethod
    def _delete_from_db(cls, db, jobid: str):
        db.execute("DELETE FROM jobs WHERE jobid=?;", (jobid, ))

    @classmethod
    def _create_uuid(cls):
        return str(uuid.uuid1())


if __name__ == "__main__":
    manager = JobManager()
    # manager.start_job("pgdb")
    for j in manager.jobs:
        print(j.to_tuple())

    print("After delete")

    manager.delete_job(jobid="f8990128-b155-11e6-ad15-dc85de763f2b")
    for j in manager.jobs:
        print(j.to_tuple())
    manager.close()


