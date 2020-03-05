from django.views.generic import TemplateView
from django.shortcuts import render
import json


class IndexView(TemplateView):
    template_name = "frontend/index.html"

class NonRelease(TemplateView):
    template_name = "frontend/non-release.html"