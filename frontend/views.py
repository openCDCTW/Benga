from django.views.generic import TemplateView


class IndexView(TemplateView):
    template_name = "frontend/index.html"


class ProfileResultView(TemplateView):
    template_name = "frontend/profile-result.html"


class ClusteringResultView(TemplateView):
    template_name = "frontend/clustering-result.html"


class TrackingResultView(TemplateView):
    template_name = "frontend/tracking-result.html"
