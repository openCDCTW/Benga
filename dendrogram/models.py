import uuid
from django.db import models


LINKAGE_CHOICES = (
    ("single", "single"),
    ("average", "average"),
)


def profiles_path(instance, filename):
    return "uploads/{0}/{1}".format(str(instance.batch_id.id), filename)


def dendrograms_png_path(instance, filename):
    return "dendrograms_png/{0}/{1}".format(str(instance.id.id),
                                            "dendrogram-" + str(instance.id.id)[0:8] + ".png")


def dendrograms_pdf_path(instance, filename):
    return "dendrograms_pdf/{0}/{1}".format(str(instance.id.id),
                                            "dendrogram-" + str(instance.id.id)[0:8] + ".pdf")


def dendrograms_svg_path(instance, filename):
    return "dendrograms_svg/{0}/{1}".format(str(instance.id.id),
                                            "dendrogram-" + str(instance.id.id)[0:8] + ".svg")


def dendrograms_emf_path(instance, filename):
    return "dendrograms_emf/{0}/{1}".format(str(instance.id.id),
                                            "dendrogram-" + str(instance.id.id)[0:8] + ".emf")


def dendrograms_newick_path(instance, filename):
    return "dendrograms_newick/{0}/{1}".format(str(instance.id.id),
                                               "dendrogram-" + str(instance.id.id)[0:8] + ".newick")


class Batch(models.Model):
    id = models.UUIDField(default=uuid.uuid4, editable=False, primary_key=True, null=False, auto_created=True)
    created = models.DateTimeField(auto_now_add=True)
    prof_num = models.SmallIntegerField(null=True)
    linkage = models.CharField(max_length=100, choices=LINKAGE_CHOICES, null=True, blank=True)


class Profile(models.Model):
    id = models.UUIDField(default=uuid.uuid4, editable=False, primary_key=True, null=False, auto_created=True)
    batch_id = models.ForeignKey(Batch, on_delete=models.CASCADE, related_name='+')
    file = models.FileField(upload_to=profiles_path, null=False)


class Dendrogram(models.Model):
    id = models.OneToOneField(Batch, on_delete=models.CASCADE, primary_key=True)
    created = models.DateTimeField(auto_now_add=True)
    linkage = models.CharField(max_length=100, choices=LINKAGE_CHOICES, null=False)
    png_file = models.FileField(upload_to=dendrograms_png_path, null=False)
    pdf_file = models.FileField(upload_to=dendrograms_pdf_path, null=False)
    svg_file = models.FileField(upload_to=dendrograms_svg_path, null=False)
    emf_file = models.FileField(upload_to=dendrograms_emf_path, null=False)
    newick_file = models.FileField(upload_to=dendrograms_newick_path, null=False)
