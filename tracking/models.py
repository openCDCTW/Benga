import uuid
from django.db import models

PROFILE_DB_CHOICES = (
    ("Vibrio_cholerae", "Vibrio cholerae"),
)


def profile_path(instance, filename):
    return "tracking/{0}/{1}".format(instance.id, "profile.tsv")


def result_path(instance, filename):
    return "tracked_results/{0}/{1}".format(str(instance.id.id), instance.json.name)


class Profile(models.Model):
    id = models.UUIDField(default=uuid.uuid4, editable=False, primary_key=True, null=False, auto_created=True)
    created = models.DateTimeField(auto_now_add=True)
    file = models.FileField(upload_to=profile_path, null=False)
    profile_db = models.CharField(max_length=100, choices=PROFILE_DB_CHOICES, null=False)


class TrackedResults(models.Model):
    id = models.OneToOneField(Profile, on_delete=models.CASCADE, primary_key=True)
    json = models.FileField(upload_to=result_path, null=False, max_length=250)
