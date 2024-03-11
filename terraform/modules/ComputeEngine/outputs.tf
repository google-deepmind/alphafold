output "connect" {
  value = "gcloud compute ssh ${google_compute_instance.my_instance.name} --project=${var.project_id} --zone=${var.zone}"
}