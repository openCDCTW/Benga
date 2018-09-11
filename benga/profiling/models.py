import uuid
from django.db import models

class UploadBatch(models.Model):
    id = models.UUIDField(default=uuid.uuid4, editable=False, primary_key=True, null=False, auto_created=True)
    created = models.DateTimeField(auto_now_add=True)


class Sequence(models.Model):
    id = models.UUIDField(default=uuid.uuid4, editable=False, primary_key=True, null=False, auto_created=True)
    batch_id = models.ForeignKey(UploadBatch, on_delete=models.CASCADE)
    filename = models.TextField(null=False)
    file = models.BinaryField(null=False, editable=True)


class Profile(models.Model):
    id = models.OneToOneField(UploadBatch, on_delete=models.CASCADE, primary_key=True)
    created = models.DateTimeField(auto_now_add=True)
    file = models.BinaryField(null=False, editable=True)
    occurrence = models.SmallIntegerField(null=False)
    database = models.TextField(null=False)


class Dendrogram(models.Model):
    id = models.OneToOneField(UploadBatch, on_delete=models.CASCADE, primary_key=True)
    created = models.DateTimeField(auto_now_add=True)
    png_file = models.BinaryField(null=False)
    pdf_file = models.BinaryField(null=False)
    svg_file = models.BinaryField(null=False)
    newick_file = models.BinaryField(null=False)