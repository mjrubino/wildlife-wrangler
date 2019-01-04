/*  Tests of how to work with occurrence data tables
Work with occs */
SELECT load_extension('mod_spatialite');

/* 1.  Copy the table  */
CREATE TABLE copy3 AS SELECT * FROM occs;

/*  2.  Add a wgs84 geometry column to the copy table*/
SELECT AddGeometryColumn('copy3', 'geom_4326', 4326, 'POINT', 'XY');

/*  3.  Add an albers geometry column to the copy table*/
ALTER TABLE 'copy3' ADD COLUMN geom_102008 BLOB;
SELECT RecoverGeometryColumn('copy3', 'geom_102008', 102008, 'POINT', 2);

/*  4.  Filter for records from june in years between 2000 and 2010 */
SELECT * FROM copy WHERE occurrenceMonth = 6;

/* Verify that the geometry is as it should be */


/*  5. Export a shapefile of june records  */
SELECT ExportSHP('copy3', 'geom_102008', 'cuckoooo', 'utf-8');
