DROP TABLE IF EXISTS `uniprot2taxid` ;
DROP TABLE IF EXISTS `uniprot_entry2accession` ;
DROP TABLE IF EXISTS `uniprot2GO` ;
DROP TABLE IF EXISTS `uniprot` ;

CREATE TABLE IF NOT EXISTS `uniprot` (
  `entry`  VARCHAR(255)           NOT NULL ,
  `source` ENUM('swissprot','trembl') NULL ,
  PRIMARY KEY (`entry`) ,
  INDEX `src_idx` (`source` ASC, `entry` ASC) )
ENGINE = InnoDB;

CREATE TABLE IF NOT EXISTS `uniprot2GO` (
  `uniprot` VARCHAR(255) NOT NULL ,
  `GO`      INT UNSIGNED NOT NULL ,
  INDEX `fk_uniprot2GO_GO2_idx` (`GO` ASC) ,
  INDEX `fk_uniprot2GO_uniprot1_idx` (`uniprot` ASC) ,
  PRIMARY KEY (`uniprot`, `GO`) ,
  CONSTRAINT `fk_uniprot2GO_GO2`
    FOREIGN KEY (`GO` )
    REFERENCES `GO` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_uniprot2GO_uniprot1`
    FOREIGN KEY (`uniprot` )
    REFERENCES `uniprot` (`entry` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;

CREATE TABLE IF NOT EXISTS `uniprot_entry2accession` (
  `entry`     VARCHAR(255) NOT NULL ,
  `accession` VARCHAR(255) NOT NULL ,
  INDEX `fk_table1_uniprot1_idx` (`entry` ASC) ,
  PRIMARY KEY (`entry`, `accession`) ,
  INDEX `access_idx` (`accession` ASC) ,
  CONSTRAINT `fk_table1_uniprot1`
    FOREIGN KEY (`entry` )
    REFERENCES `uniprot` (`entry` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;

CREATE TABLE IF NOT EXISTS `uniprot2taxid` (
  `uniprot`       VARCHAR(255) NOT NULL ,
  `taxid`         INT UNSIGNED     NULL ,
  `relationship` ENUM('is','host') NULL ,
  INDEX `fk_uniprot2taxid_uniprot1_idx` (`uniprot` ASC) ,
  INDEX `fk_uniprot2taxid_taxid1_idx` (`taxid` ASC) ,
  INDEX `relationship_idx` (`relationship` ASC, `uniprot` ASC) ,
  CONSTRAINT `fk_uniprot2taxid_uniprot1`
    FOREIGN KEY (`uniprot` )
    REFERENCES `uniprot` (`entry` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_uniprot2taxid_taxid1`
    FOREIGN KEY (`taxid` )
    REFERENCES `taxid` (`id` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;
