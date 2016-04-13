####============================================================
##  Utils for HMDB database.
##
####------------------------------------------------------------
.hmdb_column_2_xml_mapping <- c(accession="/metabolite/accession",
                                version="/metabolite/version",
                                name="/metabolite/name",
                                chem_formula="/metabolite/chemical_formula",
                                avg_molecular_weight="/metabolite/average_molecular_weight",
                                monoisotopic_molecular_weight="/metabolite/monisotopic_moleculate_weight",
                                iupac_name="/metabolite/iupac_name",
                                smiles="/metabolite/smiles",
                                inchi="/metabolite/inchi",
                                inchikey="/metabolite/inchikey")
####============================================================
##  .simpleParseHmdbXml
##
##  parse an HMDB xml file, extract some properties and return that as a data.frame.
####------------------------------------------------------------
.simpleParseHmdbXml <- function(x){
    ## Parse the XML.
    ## message("XML file", x)
    x <- xmlRoot(xmlParse(x))
    ## Check if we've got an accession!!!
    if(length(getNodeSet(x, "/metabolite/accession")) == 0)
        stop("Can not find an attribute 'accesssion' in the XML file!")
    ## From these entries I expect to get only a single value:
    basic_vals_str <- c("accession", "name", "chem_formula",
                        "iupac_name", "smiles", "inchi", "inchikey")
    strVals <- lapply(.hmdb_column_2_xml_mapping[basic_vals_str], function(z){
        res <- getNodeSet(x, z, fun=xmlValue)
        if(length(res) == 0)
            return(NA)
        return(as.character(res))
    })
    ## Process numeric entries.
    basic_vals_num <- c("version", "avg_molecular_weight", "monoisotopic_molecular_weight")
    numVals <- lapply(.hmdb_column_2_xml_mapping[basic_vals_num], function(z){
        res <- getNodeSet(x, z, fun=xmlValue)
        if(length(res) == 0)
            return(NA)
        return(as.numeric(res))
    })

    return(data.frame(strVals, numVals, stringsAsFactors=FALSE))
}

####============================================================
##  .simpleParseHmdbXmlFolder
##
##  Simply parses each xml file in the specified folder and returns
##  the extracted data as a data.frame.
####------------------------------------------------------------
.simpleParseHmdbXmlFolder <- function(x){
    xmls <- dir(x, pattern="xml", full.names=TRUE)
    message("Got ", length(xmls), " xml files; will parse all of em...", appendLF=FALSE)
    dfs <- lapply(xmls, .simpleParseHmdbXml)
    message("OK")
    res <- do.call(rbind, dfs)
    return(res)
}

####============================================================
##  .createSimpleCompoundDb
##
##  Create a simple compound identification database, here and first
##  of HMDB data.
####------------------------------------------------------------
.createSimpleCompoundDb <- function(hmdbPath=NULL, hmdbVersion=NULL, version="0.0.1",
                                    fileName="scDb.sqlite"){
    if(is.null(hmdbPath))
        stop("The path to the HMDB xml files is empty!")
    ## Metadata priming.
    met <- data.frame(name=c("Db type", "Supporting package", "Db creation date", "version",
                             "DBSCHEMAVERSION"),
                      value=c("SimpleCompoundDb", "xcmsExtensions", date(), version, "1.0"),
                      stringsAsFactors=FALSE)
    ## Initialize database:
    message("o Initialize database: ", fileName, "...", appendLF = FALSE)
    con <- dbConnect(dbDriver("SQLite"), dbname=fileName)
    message("OK")
    ## HMDB.
    if(!is.null(hmdbPath)){
        message("\n  ---- HMDB ----\n")
        xmlFiles <- dir(hmdbPath, pattern="xml", full.names=TRUE)
        if(length(xmlFiles) == 0)
            stop("No XML files found in folder: ", hmdbPath, "!")
        message("  o Processing ", length(xmlFiles), " xml files...", appendLF=FALSE)
        ## Extract the "simple" data from the xml files.
        suppressMessages(
            res <- .simpleParseHmdbXmlFolder(hmdbPath)
        )
        message("OK")
        if(nrow(res) == 0){
            stop("Could not extract any data from the xml files in folder ", hmdbPath,
                 ". Are these HMDB xml files?")
        }
        res <- cbind(res, source=rep("HMDB", nrow(res)), organism=rep("Homo sapiens", nrow(res)),
                     stringsAsFactors=FALSE)
        res <- res[, c("accession", "version", "name", "organism", "chem_formula",
                       "avg_molecular_weight", "monoisotopic_molecular_weight", "iupac_name",
                       "smiles", "inchi", "inchikey", "source")]
        res <- unique(res)
        message("  o Writing data to database...", appendLF=FALSE)
        dbWriteTable(con, name="compound_basic", res, append=TRUE)
        message("OK")
        ## Generate the metadata contents.
        met <- rbind(met, c("db:HMDB", hmdbVersion))
    }
    ## Metadata
    message("o Processing metadata...", appendLF=FALSE)
    dbWriteTable(con, "metadata", met, append=FALSE, overwrite=TRUE)
    message("OK")
    ## Indices:
    message("o Creating indices...", appendLF=FALSE)
    dbGetQuery(con, "create index acc_idx on compound_basic (accession);")
    dbGetQuery(con, "create index avg_weight_idx on compound_basic (avg_molecular_weight);")
    dbGetQuery(con, "create index mono_weight_idx on compound_basic (monoisotopic_molecular_weight);")
    message("OK")
    dbDisconnect(con)
    return(fileName)
}
