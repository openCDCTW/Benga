import hashlib
import uuid
from django.db import models


def sequences_path(instance, filename):
    filename = hashlib.sha256(instance.file.read()).hexdigest()
    return "sequences/{0}/{1}".format(instance.batch_id.id, filename + ".fa")


def profiles_path(instance, filename):
    return "profiles/{0}/{1}".format(instance.id, "profile-" + str(instance.id)[0:8] + ".tsv")


def dendrograms_png_path(instance, filename):
    return "dendrograms_png/{0}/{1}".format(instance.id, "dendrogram-" + str(instance.id)[0:8] + ".png")


def dendrograms_pdf_path(instance, filename):
    return "dendrograms_pdf/{0}/{1}".format(instance.id, "dendrogram-" + str(instance.id)[0:8] + ".pdf")


def dendrograms_svg_path(instance, filename):
    return "dendrograms_svg/{0}/{1}".format(instance.id, "dendrogram-" + str(instance.id)[0:8] + ".svg")


def dendrograms_newick_path(instance, filename):
    return "dendrograms_newick/{0}/{1}".format(instance.id, "dendrogram-" + str(instance.id)[0:8] + ".newick")


class UploadBatch(models.Model):
    id = models.UUIDField(default=uuid.uuid4, editable=False, primary_key=True, null=False, auto_created=True)
    created = models.DateTimeField(auto_now_add=True)


class Sequence(models.Model):
    id = models.UUIDField(default=uuid.uuid4, editable=False, primary_key=True, null=False, auto_created=True)
    batch_id = models.ForeignKey(UploadBatch, on_delete=models.CASCADE)
    file = models.FileField(upload_to=sequences_path, null=False)


class Profile(models.Model):
    id = models.OneToOneField(UploadBatch, on_delete=models.CASCADE, primary_key=True)
    created = models.DateTimeField(auto_now_add=True)
    file = models.FileField(upload_to=profiles_path, null=False)
    occurrence = models.SmallIntegerField(null=False)
    database = models.TextField(null=False)


class Dendrogram(models.Model):
    id = models.OneToOneField(UploadBatch, on_delete=models.CASCADE, primary_key=True)
    created = models.DateTimeField(auto_now_add=True)
    png_file = models.FileField(upload_to=dendrograms_png_path, null=False)
    pdf_file = models.FileField(upload_to=dendrograms_pdf_path, null=False)
    svg_file = models.FileField(upload_to=dendrograms_svg_path, null=False)
    newick_file = models.FileField(upload_to=dendrograms_newick_path, null=False)