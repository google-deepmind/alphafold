locals {
  service_account_email = var.service_account_email == null ? data.google_compute_default_service_account.default[0].email : var.service_account_email
}

data "google_compute_default_service_account" "default" {
  count = var.service_account_email==null?1:0
}

resource "google_compute_instance" "my_instance" {
  name         = var.instance_name
  project = var.project_id 
  machine_type = var.machine_type
  zone         = var.zone

  boot_disk {
    initialize_params {
      image = var.image
      labels = var.labels
    }
  }
   guest_accelerator {
    type = var.accelerator_type
    count = var.accelerator_count
   }
    network_interface {
        network = var.network_name
        subnetwork = var.subnet_name 
    }

    service_account {
    # Google recommends custom service accounts that have cloud-platform scope and permissions granted via IAM Roles.
    email  = local.service_account_email
    scopes = ["cloud-platform"]
   }

    scheduling {
    preemptible        = false
    automatic_restart  = true
    on_host_maintenance = "TERMINATE"  # Disable live migration
  }
    metadata_startup_script = var.startup_script

  }
