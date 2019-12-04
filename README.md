## Purpose
The abundance of wildlife occurrence datasets that are currently accessible can be valuable for efforts such as species distribution modeling and range delineation.  However, the task of downloading and cleaning occurrence records is more complex than it may seem at first consideration given errors and uncertainties in data.  This repository is a draft framework for collecting and filtering occurrence data that is widely available through API's.  

## Framework
Data is requested from occurrence dataset API's and filtered according to species- and request-specific parameters.  Filtered occurrence records are saved in a database.  The details of species-concepts and filter parameter sets are stored in a database for use and reference.  Additionally, jupyter notebooks are created that describe the filtered datasets for the sake of documentation and for decision making about filter parameterization.

## Features
This framework is designed to meet several important criteria in order to provide summaries that can be interpreted at face value with high confidence or trusted when feeding into analyses and evaluations.
* Automation -- the volume of data and species involved necessitates the processes be automated. Automation also reduces subjectivity in decision making, enables thorough documentation, and ensures repeatability.

* High confidence -- criteria and filters can be set in this framework that produce high confidence in results.

* Detailed parameterization -- data requests and record filters can be parameterized on a per-species and per-event basis. Rules do not have to be applied generally across large numbers of species or evaluations.

* Open source -- processes are coded in Python 3 and SQL and use sqlite3, which comes with Python 3, for spatial queries.

* Transparency through documentation -- summaries and models using empirical data still involve subjectivity in the form of choices regarding parameters and rules for handling/filtering/cleaning data.  This framework is meant to document those choices.  It identifies several filter criteria that are especially relevant to species-occurrence records.
  * Occurrence year -- species' distributions and taxonomy can change over time.
  * Occurrence month -- especially relevant regarding migratory species.
  * Coordinate uncertainty and issues -- this varies immensely among records and some records have known issues that limit their value.
  * Detection distance -- different species may be sampled differently in ways that introduce location uncertainty.  Small mammals that are captured in traps can confidently be assigned to the trap's location, but loud-singing birds detected in an auditory survey could be hundreds of meters away from the observer and thus, the coordinate associated with the record.  
  * Occurrence records are circles -- although records are recorded as x,y coordinates, the coordinate uncertainty and detection distance issues described above require that they be treated as circles with centers at the x,y coordinate and radius equal to the detection distance plus the coordinate uncertainty.  Buffering the points accounts for this.  
  * Occurrence data sets are dynamic -- some datasets enable data contributors to go back and edit attributes of occurrence records.  In addition, historic records may be added that change the set of records associated with a past time period.  That is to say that a query of years past run today may be different than the same query run tomorrow.  This represents a challenge for provenance.  The method I employed to handle this is to document data request parameters and post-request filter sets as uniquely identifiable objects that are stored and documented in a database ('parameters.sqlite').  Records that pass through filters are also stored with geometry in a separate sqlite database.  This could be moved to sql server or postgresql later.
  * Species concepts are dynamic -- taxonomic classifications are constantly being scrutinized and are revised annually.  As a result, different projects may have used different concepts for the same species name (common or scientific).  In many cases, this is not a big problem because the species is unique and easily identifiable and taxonomic changes regard names only.  More problematic are cases where genetic studies have identified species that are nearly identical physiologically and were once considered a single species.  Such cases often reveal geographic patterns in species occurrence that may be used by some observers as a basis for identification of individuals in the field thus introducing circularity. For example, a bird watcher decides an individual's identity in part based on which species is supposed to occur in the area *or* a species' range is revised and eBird changes records of the old species concept from the area where it is now known not to occur to the correct species.  Dynamic species-concepts are a daunting challenge, but can hopefully be handled with scrutinizing taxonomic crosswalks.  In this framework, species-concepts are documented in 'parameters.sqlite' with columns for the years the concept was valid and a geometry column where a polygon of potential occurrence could be recorded.  This whole topic needs more work, but the BCB Taxa Information Registry can contribute a lot.
  * Data sensivity regarding poaching risk -- More a challenge than criteria, is the issue of occurrence records for many species, especially herpetofauna, being sensitive data due to the threat of poaching.  Individuals of rare and imperiled species are regularly captured in the wild for sale on black markets, sometimes providing hundreds of dollars per individual.  Poachers are known to review and interpret scientific data in order to determine the locations of populations that they can collect from.  This creates a serious challenge for those managing and summarizing locational data.  Managers of some databases "fuzz" or buffer points to coarsen the information to a level that isn't useful to poachers, which also limits the usefulness of data for conservation assessments.  Not only does this issues create limitations on the user's end, it also creates restrictions on how data providers can/should share and store records.  Data in the hands of federal agencies is subject to Freedom of Information Act requests.  This issue may shape aspects of the framework.  For example, it may be necessary to make this technology deployable to individual users and include a portal for data sets that are available to the user but remain unavailable from public facing or federal data sets.  

## Recent changes
* Added ability to limit requests to within geometries.
* Handling duplicate records on the bases of x y coordinats and date-time
* Added ability to specify a limiting extent for a species concept.
* Increased process speed by 2,300%.

## Coming soon
* Handling requests of > 200,000 records.
* Ability to incorporate records directly from a copy of the eBird data set that user has downloaded.
* Generalization of framework to better facilitate multiple users.  Parameters.sqlite will be replaced with a template that user can build from.  Each user will access a local copy of parameters.sqlite filled out to for their needs.

## Inputs
Data is gathered from catalogs and databases through API's so there are few inputs.  However, the 'parameters.sqlite' database is needed, which includes tables for species-concepts, data request parameters, and post-request filtering criteria.

GBIF is currently the only dataset currently used but others can/will be added later including eBird.

## Outputs
On a per-species basis
* A database of filtered species-occurrence records.  The format supports display in QGIS and other GIS.
* Notebooks that describe decisions made and the data acquired.

## Constraints
* Queries returning > 200,000 records will fail.

## Dependencies
Python 3 and numerous packages including sqlite3 with the spatialite extension.  An environment can be created from the ENVIRONMENT.yml file included in this repository.

## Code
All code is included in this repository.  Runtimes of discrete tasks made grouping code into separate functions preferable.  
