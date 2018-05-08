DROP TABLE IF EXISTS `GO_extended` ;
DROP TABLE IF EXISTS `GO_alternative` ;
DROP TABLE IF EXISTS `GO_relationship` ;
DROP TABLE IF EXISTS `GO` ;

CREATE TABLE IF NOT EXISTS `GO` (
  `nid`       INT UNSIGNED  NOT NULL AUTO_INCREMENT ,
  `id`        VARCHAR(255)  NOT NULL ,
  `name`      VARCHAR(510)      NULL ,
  `namespace` ENUM('B','M','C') NULL COMMENT 'B: biological_process\nM: molecular_function\nC: cellular_component' ,
  `obsolete`  TINYINT(1)        NULL DEFAULT 0 ,
  UNIQUE INDEX `id_UNIQUE` (`id` ASC) ,
  INDEX `namespace_idx` (`namespace` ASC, `id` ASC) ,
  PRIMARY KEY (`nid`) ,
  INDEX `go_id` (`nid` ASC) ,
  INDEX `name_idx` (`namespace` ASC, `name` ASC) )
ENGINE = InnoDB;

CREATE TABLE IF NOT EXISTS `GO_extended` (
  `parent_nid` INT UNSIGNED NOT NULL ,
  `child_nid`  INT UNSIGNED NOT NULL ,
  INDEX `fk_GO_extended_GO2_idx` (`child_nid` ASC) ,
  INDEX `parent_idx` (`parent_nid` ASC) ,
  PRIMARY KEY (`parent_nid`, `child_nid`) ,
  CONSTRAINT `fk_GO_extended_GO1`
    FOREIGN KEY (`parent_nid` )
    REFERENCES `GO` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_GO_extended_GO2`
    FOREIGN KEY (`child_nid` )
    REFERENCES `GO` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;

CREATE TABLE IF NOT EXISTS `GO_alternative` (
  `GO_nid`      INT UNSIGNED NOT NULL ,
  `alternative` VARCHAR(255) NOT NULL ,
  INDEX `fk_GO_alternative_GO1_idx` (`GO_nid` ASC) ,
  INDEX `alternative_idx` (`alternative` ASC) ,
  CONSTRAINT `fk_GO_alternative_GO1`
    FOREIGN KEY (`GO_nid` )
    REFERENCES `GO` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;

CREATE TABLE IF NOT EXISTS `GO_relationship` (
  `source_nid`      INT UNSIGNED NOT NULL ,
  `destination_nid` INT UNSIGNED NOT NULL ,
  `relationship`    ENUM('part_of','regulates','negatively_regulates','positively_regulates','has_part','occurs_in','results_in') NOT NULL ,
  INDEX `fk_GO_relationship_GO1_idx` (`source_nid` ASC) ,
  INDEX `fk_GO_relationship_GO2_idx` (`destination_nid` ASC) ,
  PRIMARY KEY (`source_nid`, `destination_nid`) ,
  INDEX `relationship_idx` (`relationship` ASC, `source_nid` ASC, `destination_nid` ASC) ,
  CONSTRAINT `fk_GO_relationship_GO1`
    FOREIGN KEY (`source_nid` )
    REFERENCES `GO` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_GO_relationship_GO2`
    FOREIGN KEY (`destination_nid` )
    REFERENCES `GO` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;