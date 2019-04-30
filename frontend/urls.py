from django.urls import path
from . import views

urlpatterns = [
    path('', views.IndexView.as_view()),
    path('profile/<uuid:pk>/', views.ProfileResultView.as_view()),
    path('clustering-results/<uuid:pk>/', views.ClusteringResultView.as_view()),
    path('tracking-results/<uuid:pk>/', views.TrackingResultView.as_view()),
]
