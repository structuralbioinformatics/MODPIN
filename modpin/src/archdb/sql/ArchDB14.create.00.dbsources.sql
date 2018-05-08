/*
*
SQL CREATOR FOR ARCHDB 2014 (ArchDB14)
SECTION 00. SOURCE DATABASES
jbonet @ oliva's lab
*
*/

DROP   TABLE IF     EXISTS `source_databases` ;
CREATE TABLE IF NOT EXISTS `source_databases` (
    `dbname`   VARCHAR(45) NOT NULL ,
    `dbsource` VARCHAR(255)    NULL ,
    `date`     DATE            NULL ,
    PRIMARY KEY (`dbname`) )
ENGINE = InnoDB;

DROP function IF EXISTS `strSplit`;

DELIMITER $$
CREATE FUNCTION strSplit(x varchar(255), delim varchar(12), pos int) returns varchar(255)
return replace(substring(substring_index(x, delim, pos), length(substring_index(x, delim, pos - 1)) + 1), delim, '');$$

DELIMITER ;