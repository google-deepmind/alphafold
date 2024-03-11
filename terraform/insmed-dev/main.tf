module "computeEngine" {
  source            = "../modules/ComputeEngine"
  project_id      = var.project_id
  count             = length(var.ComputeEngine)
  instance_name     = var.ComputeEngine[count.index].instance_name
  machine_type      = var.ComputeEngine[count.index].machine_type
  zone              = var.ComputeEngine[count.index].zone
  network_name      = var.ComputeEngine[count.index].network_name
  subnet_name       = var.ComputeEngine[count.index].subnet_name
  startup_script    = var.ComputeEngine[count.index].startup_script
  labels            = var.ComputeEngine[count.index].labels
  image             = var.ComputeEngine[count.index].image
  accelerator_type  = var.ComputeEngine[count.index].accelerator_type
  accelerator_count = var.ComputeEngine[count.index].accelerator_count
  service_account_email = var.ComputeEngine[count.index].service_account_email
}

