# Database functions for HGIS


# Create SQLite database for HGIS



create_hgis_db <- function() {
  # Create database
  hgisdb <- DBI::dbConnect(RSQLite::SQLite(), here::here("data/hgis.sqlite"))
  
  # Set up tables
  schema <- c("CREATE TABLE hgis_raw (
               runtime TEXT NOT NULL,
               sample_id INTEGER,
               wheel_pos_id INTEGER,
               ok_calc INTEGER NOT NULL,
               mst_num INTEGER NOT NULL,
               le12c REAL,
               le13c REAL,
               he12c REAL,
               he13c REAL,
               cnt14c INTEGER,
               he13_12 REAL,
               corr_14_12 REAL,
               sig_14_12 REAL,
               ltcorr REAL,
               norm_ratio REAL,
               sig_norm_ratio REAL,
               norm_del13c REAL,
               sig_norm_del13c REAL
             );",
              
             "CREATE TABLE hgis_samples (
               sample_id INTEGER PRIMARY KEY,
               rec_num INTEGER,
               name TEXT,
               carbonate_mass REAL
               helium_added REAL
             );",
              
             "CREATE TABLE hgis_results (
               result_id INTEGER PRIMARY KEY,
               sample_id INTEGER,
               runtime TEXT,
               sample_type TEXT,
               num_runs INTEGER,
               tot_runs INTEGER,
               norm_ratio REAL,
               int_err REAL,
               ext_err REAL,
               del_13c REAL,
               sig_del_13c REAL,
               fm_corr REAL,
               sig_fm_corr REAL,
               lg_blk_fm,
               sig_lg_blk_fm
             );",
             
             "CREATE TABLE wheel_pos (
               wheel_pos_id INTEGER PRIMARY KEY,
               sample_id INTEGER,
               wheel_id TEXT,
               pos INTEGER
             );",
             
             "CREATE TABLE hgis_conditions (
               hgis_conditions_id INTEGER PRIMARY KEY,
               sample_id INTEGER,
               initial_current REAL,
               base_current REAL,
               source_pressure REAL,
               helium_flow REAL,
               cap_length REAL,
               cap_id REAL
             );")
  
  purrr::walk(schema, function(x) DBI::dbExecute(hgisdb, x))
}


#' Insert a table of samples into HGIS DB
#'
#' @param data A dataframe of sample info in `hgis_samples` format.
#' @param con A DBI connection object for the HGIS DB
#'
#' @return
#' @export
#'
insert_samples <- function(data, con = hgisdb) {
  
  next_sample_id <- dbGetQuery(con, "SELECT MAX(sample_id) FROM hgis_samples") %>% 
    pull(sample_id) + 1
  
  data$sample_id <- next_sample_id:(nrow(data) + next_sample_id)
  dbAppendTable(con, "hgis_samples", data)
}

# Insert raw data

format_insert <- function(data) {
  data <- data %>% 
    mutate(runtime = as.character(ts),
           ok_calc = !as.numeric(outlier),
           ltcorr = Cycles/10,
           sig_norm_ratio = ce/normFm,
           le12c = le12C * 1E-6,
           he12c = he12C * 1E-6) %>% 
    select(ok_calc, runtime, mst_num = Meas,
           le12c, le13c = le13C, he12c,
           he13c = he13C, cnt14c = CntTotGT, he13_12 = X13.12he,
           corr_14_12 = cor1412he, sig_14_12 = ce, ltcorr,
           norm_ratio = normFm, sig_norm_ratio)
}

insert_raw <- function(data, con = hgisdb) {
  dbAppendTable(con, "hgis_raw", data)
}

# Insert sample data


# Insert results

# Sample conditions?


# Get wheel data

# Get secondary data

