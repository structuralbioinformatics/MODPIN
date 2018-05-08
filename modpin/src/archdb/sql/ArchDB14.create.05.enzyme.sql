DROP TABLE IF EXISTS `uniprot2enzyme` ;
DROP TABLE IF EXISTS `enzyme_compound` ;
DROP TABLE IF EXISTS `enzyme_cofactor` ;
DROP TABLE IF EXISTS `enzyme_reaction` ;
DROP TABLE IF EXISTS `enzyme_reaction_compounds` ;
DROP TABLE IF EXISTS `enzyme_extended` ;
DROP TABLE IF EXISTS `enzyme_transfered` ;
DROP TABLE IF EXISTS `enzyme` ;

CREATE TABLE IF NOT EXISTS `enzyme` (
  `id`          VARCHAR(45)  NOT NULL ,
  `description` VARCHAR(255) NOT NULL ,
  `level`       TINYINT UNSIGNED NULL DEFAULT NULL ,
  `parent`      VARCHAR(45)      NULL DEFAULT NULL ,
  PRIMARY KEY (`id`) ,
  INDEX `fk_enzyme_enzyme_idx` (`parent` ASC) ,
  CONSTRAINT `fk_enzyme_enzyme`
    FOREIGN KEY (`parent` )
    REFERENCES `enzyme` (`id` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;

CREATE TABLE IF NOT EXISTS `enzyme_compound` (
  `nid`  INT UNSIGNED                NOT NULL AUTO_INCREMENT ,
  `name` VARCHAR(510)                NOT NULL ,
  `type` ENUM('cofactor','compound') NOT NULL ,
  PRIMARY KEY (`nid`) ,
  UNIQUE INDEX `name_UNIQUE` (`name` ASC, `type` ASC) )
ENGINE = InnoDB;

CREATE TABLE IF NOT EXISTS `enzyme_cofactor` (
  `ec`       VARCHAR(45)  NOT NULL ,
  `cofactor` INT UNSIGNED NOT NULL ,
  INDEX `fk_enzyme_cofactor_enzyme_idx` (`ec` ASC) ,
  PRIMARY KEY (`ec`, `cofactor`) ,
  INDEX `cofactor_index` (`cofactor` ASC) ,
  CONSTRAINT `fk_enzyme_cofactor_enzyme`
    FOREIGN KEY (`ec` )
    REFERENCES `enzyme` (`id` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_enzyme_cofactor_enzyme_compound1`
    FOREIGN KEY (`cofactor` )
    REFERENCES `enzyme_compound` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;

CREATE TABLE IF NOT EXISTS `enzyme_reaction` (
  `nid` INT UNSIGNED NOT NULL AUTO_INCREMENT ,
  `ec`  VARCHAR(45)  NOT NULL ,
  `description` VARCHAR(510) NOT NULL ,
  INDEX `fk_enzyme_reaction_enzyme_idx` (`ec` ASC) ,
  INDEX `reaction_index` (`nid` ASC) ,
  CONSTRAINT `fk_enzyme_reaction_enzyme`
    FOREIGN KEY (`ec` )
    REFERENCES `enzyme` (`id` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;

CREATE TABLE IF NOT EXISTS `enzyme_reaction_compounds` (
  `reaction`    INT UNSIGNED                NOT NULL ,
  `compound`    INT UNSIGNED                NOT NULL ,
  `type`        ENUM('substrate','product') NOT NULL ,
  `cardinality` TINYINT UNSIGNED            NOT NULL DEFAULT 1 ,
  INDEX `fk_reaction_compounds_enzyme_reaction_idx` (`reaction` ASC) ,
  INDEX `compound_idx` (`compound` ASC, `cardinality` ASC) ,
  INDEX `type_comp_idx` (`type` ASC, `compound` ASC, `cardinality` ASC) ,
  CONSTRAINT `fk_reaction_compounds_enzyme_reaction`
    FOREIGN KEY (`reaction` )
    REFERENCES `enzyme_reaction` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_enzyme_reaction_compounds_enzyme_compound1`
    FOREIGN KEY (`compound` )
    REFERENCES `enzyme_compound` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;

CREATE TABLE IF NOT EXISTS `enzyme_extended` (
  `parent_id` VARCHAR(45) NOT NULL ,
  `child_id`  VARCHAR(45) NOT NULL ,
  INDEX `fk_enzyme_extended_enzyme1_idx` (`parent_id` ASC) ,
  INDEX `fk_enzyme_extended_enzyme2_idx` (`child_id` ASC) ,
  CONSTRAINT `fk_enzyme_extended_enzyme1`
    FOREIGN KEY (`parent_id` )
    REFERENCES `enzyme` (`id` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_enzyme_extended_enzyme2`
    FOREIGN KEY (`child_id` )
    REFERENCES `enzyme` (`id` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;

CREATE TABLE IF NOT EXISTS `enzyme_transfered` (
  `original_id`   VARCHAR(45) NOT NULL ,
  `transfered_id` VARCHAR(45) NOT NULL ,
  INDEX `fk_enzyme_transfered_enzyme1_idx` (`original_id` ASC) ,
  INDEX `fk_enzyme_transfered_enzyme2_idx` (`transfered_id` ASC) ,
  CONSTRAINT `fk_enzyme_transfered_enzyme1`
    FOREIGN KEY (`original_id` )
    REFERENCES `enzyme` (`id` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_enzyme_transfered_enzyme2`
    FOREIGN KEY (`transfered_id` )
    REFERENCES `enzyme` (`id` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;

CREATE TABLE IF NOT EXISTS `uniprot2enzyme` (
  `uniprot` VARCHAR(255) NOT NULL ,
  `enzyme`  VARCHAR(45)  NOT NULL ,
  INDEX `fk_uniprot2enzyme_uniprot1_idx` (`uniprot` ASC) ,
  INDEX `fk_uniprot2enzyme_enzyme1_idx` (`enzyme` ASC) ,
  PRIMARY KEY (`uniprot`, `enzyme`) ,
  CONSTRAINT `fk_uniprot2enzyme_uniprot1`
    FOREIGN KEY (`uniprot` )
    REFERENCES `uniprot` (`entry` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_uniprot2enzyme_enzyme1`
    FOREIGN KEY (`enzyme` )
    REFERENCES `enzyme` (`id` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;
