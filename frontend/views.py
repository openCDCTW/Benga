from django.views.generic import TemplateView
from django.shortcuts import render
import json


class IndexView(TemplateView):
    template_name = "frontend/index.html"


class ProfileResultView(TemplateView):
    template_name = "frontend/profile-result.html"

    def get(self, request, jobid):
        return render(request, self.template_name, {"jobid": json.dumps(jobid)})


class ClusteringResultView(TemplateView):
    template_name = "frontend/clustering-result.html"

    def get(self, request, jobid):
        return render(request, self.template_name, {"jobid": json.dumps(jobid)})


class TrackingResultView(TemplateView):
    template_name = "frontend/tracking-result.html"

    def get(self, request, jobid):
        return render(request, self.template_name, {"jobid": json.dumps(jobid)})
