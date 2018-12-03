import uuid
from django.db import models


def sequences_path(instance, filename):
    return "uploaded_sequences/{0}/{1}".format(instance.batch_id.id, instance.file.name)


def profiles_path(instance, filename):
    return "profiles/{0}/{1}".format(str(instance.id.id),
                                     "profiles-" + str(instance.id.id)[0:8] + ".tsv")


def zip_path(instance, filename):
    return "profiles/{0}/{1}".format(str(instance.id.id),
                                     "individual-profiles-" + str(instance.id.id)[0:8] + ".zip")


class Batch(models.Model):
    id = models.UUIDField(default=uuid.uuid4, editable=False, primary_key=True, null=False, auto_created=True)
    created = models.DateTimeField(auto_now_add=True)


class Sequence(models.Model):
    id = models.UUIDField(default=uuid.uuid4, editable=False, primary_key=True, null=False, auto_created=True)
    batch_id = models.ForeignKey(Batch, on_delete=models.CASCADE)
    file = models.FileField(upload_to=sequences_path, null=False)


class Profile(models.Model):
    id = models.OneToOneField(Batch, on_delete=models.CASCADE, primary_key=True)
    created = models.DateTimeField(auto_now_add=True)
    file = models.FileField(upload_to=profiles_path, null=False)
    zip = models.FileField(upload_to=zip_path, null=False)
    occurrence = models.SmallIntegerField(null=False)
    database = models.TextField(null=False)
