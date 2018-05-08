DROP TABLE IF EXISTS `taxid_old` ;
DROP TABLE IF EXISTS `taxid` ;

CREATE TABLE IF NOT EXISTS `taxid` (
  `id`     INT UNSIGNED NOT NULL ,
  `name`   VARCHAR(255) NOT NULL COMMENT 'scientific name' ,
  `parent` INT UNSIGNED NOT NULL ,
  `rank`   ENUM('no rank','superkingdom','kingdom','subkingdom','superphylum','phylum','subphylum','superclass','class','subclass','infraclass','superorder','order','suborder','infraorder','parvorder','superfamily','family','subfamily','tribe','subtribe','genus','subgenus','species group','species subgroup','species','subspecies','varietas','forma') NOT NULL ,
  PRIMARY KEY (`id`) ,
  UNIQUE INDEX `id_UNIQUE` (`id` ASC) ,
  INDEX `parent_idx` (`parent` ASC, `id` ASC) ,
  INDEX `rank_idx` (`rank` ASC, `id` ASC, `parent` ASC) )
ENGINE = InnoDB;

CREATE TABLE IF NOT EXISTS `taxid_old` (
  `oldid` INT UNSIGNED NOT NULL ,
  `taxid` INT UNSIGNED     NULL DEFAULT NULL ,
  PRIMARY KEY (`oldid`) ,
  INDEX `fk_oldtaxid_taxid1_idx` (`taxid` ASC) ,
  CONSTRAINT `fk_oldtaxid_taxid1`
    FOREIGN KEY (`taxid` )
    REFERENCES `taxid` (`id` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;