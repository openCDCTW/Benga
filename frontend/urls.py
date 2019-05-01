from django.urls import path
from . import views

urlpatterns = [
    path('', views.IndexView.as_view()),
    path('profile/<uuid:jobid>/', views.ProfileResultView.as_view()),
    path('clustering-results/<uuid:jobid>/', views.ClusteringResultView.as_view()),
    path('tracking-results/<uuid:jobid>/', views.TrackingResultView.as_view()),
]
