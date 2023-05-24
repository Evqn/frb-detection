
from django.urls import path
from . import views


urlpatterns = [
    path('', views.redirect_view), 
    path('home', views.data), 
    path('home/scan<int:scan_num>/', views.scan),
    path('home/scan<int:scan_num>/image<int:image_num>/', views.image)
]
