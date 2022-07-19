# These are some basic wrappers to provide support for keyring access and
# control for the project.
accessKey <- function(service) {

  existsQ <- tryCatch(
    {
      key <- keyring::key_get(service)
    },
    error = function(cond) {
      message(paste("The key, ",
                 service,
                    ", does not appear to exist. You will need to create ",
                    "the key with createKey(username, password)."))
    },
    finally = {
      # key <- keyring::key_list(service)
      # key$password <- keyring::key_get(service)
      return(key)
    }
  )

  keyring::default_backend()
  if ( !keyring::has_keyring_support() ) {
    stop("The operating system does not support keyring access.")
  }
  if ( !keyring::keyring_is_locked(keyring = NULL) ) {
    keyring::keyring_unlock(keyring = NULL, password = NULL)
  }
  return(existsQ)
}

createKey <- function(service, username = NULL) {
  # createKey will create a new key with the specified service name and username
  # in the default OS keychain..
  tryCatch(
    {
      keyring::key_get(service)
    },
    error = function(cond) {
      keyring::key_set(service, username)
    }
  )
}

updateKey <- function(service, username = NULL) {
  tryCatch(
    {
      keyring::key_get(service)
    },
    error = function(cond) {
      message(paste("The service, ",service,", does not exists."))
    },
    finally = {
      keyring::key_set(service, username)
    }
  )
}
