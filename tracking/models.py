import uuid
from django.db import models


def sequences_path(instance, filename):
    return "tracking/{0}/{1}".format(instance.id, instance.file.name)


def result_path(instance, filename):
    return "tracked_results/{0}/{1}".format(instance.id, instance.json.name)


class Sequence(models.Model):
    id = models.UUIDField(default=uuid.uuid4, editable=False, primary_key=True, null=False, auto_created=True)
    created = models.DateTimeField(auto_now_add=True)
    file = models.FileField(upload_to=sequences_path, null=False)


class TrackedResults(models.Model):
    id = models.OneToOneField(Sequence, on_delete=models.CASCADE, primary_key=True)
    json = models.FileField(upload_to=result_path, null=False)


class Tracking(models.Model):
    ALLELE_DB_CHOICES = (
        ("Vibrio_cholerae", "Vibrio cholerae"),
    )
    PROFILE_DB_CHOICES = (
        ("vibrio-profiles", "Vibrio cholerae"),
    )
    id = models.OneToOneField(Sequence, on_delete=models.CASCADE, primary_key=True)
    allele_db = models.CharField(max_length=100, choices=ALLELE_DB_CHOICES, null=False)
    profile_db = models.CharField(max_length=100, choices=PROFILE_DB_CHOICES, null=False)
