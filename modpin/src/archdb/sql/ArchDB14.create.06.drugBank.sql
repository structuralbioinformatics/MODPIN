DROP TABLE IF EXISTS `drugBank_target_action` ;
DROP TABLE IF EXISTS `drugBank_target` ;
DROP TABLE IF EXISTS `drugBank_synonym` ;
DROP TABLE IF EXISTS `drugBank_secondary_accession_number` ;
DROP TABLE IF EXISTS `drugBank_group` ;
DROP TABLE IF EXISTS `drugBank_external_identifiers` ;
DROP TABLE IF EXISTS `drugBank_category` ;
DROP TABLE IF EXISTS `drugBank` ;

CREATE  TABLE IF NOT EXISTS `drugBank` (
  `drugbank_id` VARCHAR(255) NOT NULL ,
  `name` VARCHAR(255) NOT NULL ,
  `type` VARCHAR(45) NOT NULL ,
  `cas_number` VARCHAR(255) NULL DEFAULT NULL ,
  `description` TEXT NULL DEFAULT NULL ,
  `hydrophobicity` FLOAT NULL DEFAULT NULL ,
  `isoelectric_point` FLOAT NULL DEFAULT NULL ,
  `molecular_weight` VARCHAR(45) NULL DEFAULT NULL ,
  `molecular_formula` VARCHAR(255) NULL DEFAULT NULL ,
  PRIMARY KEY (`drugbank_id`) )
ENGINE = InnoDB;

CREATE  TABLE IF NOT EXISTS `drugBank_category` (
  `drugbank_id` VARCHAR(255) NOT NULL ,
  `drugbank_category` VARCHAR(45) NOT NULL ,
  INDEX `fk_drugBank_category_drugBank1_idx` (`drugbank_id` ASC) ,
  CONSTRAINT `fk_drugBank_category_drugBank1`
    FOREIGN KEY (`drugbank_id` )
    REFERENCES `drugBank` (`drugbank_id` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;

CREATE  TABLE IF NOT EXISTS `drugBank_external_identifiers` (
  `drugbank_id` VARCHAR(255) NOT NULL ,
  `identifier` VARCHAR(255) NOT NULL ,
  `resource` VARCHAR(255) NOT NULL ,
  INDEX `fk_drugBank_external_identifiers_drugBank1_idx` (`drugbank_id` ASC) ,
  CONSTRAINT `fk_drugBank_external_identifiers_drugBank1`
    FOREIGN KEY (`drugbank_id` )
    REFERENCES `drugBank` (`drugbank_id` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;

CREATE  TABLE IF NOT EXISTS `drugBank_group` (
  `drugbank_id` VARCHAR(255) NOT NULL ,
  `drugbank_group` VARCHAR(45) NOT NULL ,
  INDEX `fk_drugBank_group_drugBank1_idx` (`drugbank_id` ASC) ,
  CONSTRAINT `fk_drugBank_group_drugBank1`
    FOREIGN KEY (`drugbank_id` )
    REFERENCES `drugBank` (`drugbank_id` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;

CREATE  TABLE IF NOT EXISTS `drugBank_secondary_accession_number` (
  `drugbank_id` VARCHAR(255) NOT NULL ,
  `secondary_accession_number` VARCHAR(45) NOT NULL ,
  INDEX `fk_drugBank_secondary_accession_number_drugBank1_idx` (`drugbank_id` ASC) ,
  CONSTRAINT `fk_drugBank_secondary_accession_number_drugBank1`
    FOREIGN KEY (`drugbank_id` )
    REFERENCES `drugBank` (`drugbank_id` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;

CREATE  TABLE IF NOT EXISTS `drugBank_synonym` (
  `drugbank_id` VARCHAR(255) NOT NULL ,
  `drugbank_synonym` VARCHAR(45) NOT NULL ,
  INDEX `fk_drugBank_synonym_drugBank1_idx` (`drugbank_id` ASC) ,
  CONSTRAINT `fk_drugBank_synonym_drugBank1`
    FOREIGN KEY (`drugbank_id` )
    REFERENCES `drugBank` (`drugbank_id` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;

CREATE  TABLE IF NOT EXISTS `drugBank_target` (
  `drugbank_id` VARCHAR(255) NOT NULL ,
  `uniprot` VARCHAR(255) NOT NULL ,
  `known_action` ENUM('yes','no','unknown') NULL ,
  `target_type` ENUM('target','enzyme','transporter','carrier') NULL ,
  INDEX `fk_drugBank_target_drugBank1_idx` (`drugbank_id` ASC) ,
  INDEX `fk_drugBank_target_uniprot1_idx` (`uniprot` ASC) ,
  CONSTRAINT `fk_drugBank_target_drugBank1`
    FOREIGN KEY (`drugbank_id` )
    REFERENCES `drugBank` (`drugbank_id` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_drugBank_target_uniprot1`
    FOREIGN KEY (`uniprot` )
    REFERENCES `uniprot` (`entry` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;

CREATE  TABLE IF NOT EXISTS `drugBank_target_action` (
  `drugbank_id` VARCHAR(255) NOT NULL ,
  `uniprot` VARCHAR(255) NOT NULL ,
  `action` VARCHAR(255) NULL ,
  INDEX `fk_drugBank_target_action_drugBank1_idx` (`drugbank_id` ASC) ,
  INDEX `fk_drugBank_target_action_uniprot1_idx` (`uniprot` ASC) ,
  INDEX `action_idx` (`action` ASC, `drugbank_id` ASC) ,
  CONSTRAINT `fk_drugBank_target_action_drugBank1`
    FOREIGN KEY (`drugbank_id` )
    REFERENCES `drugBank` (`drugbank_id` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_drugBank_target_action_uniprot1`
    FOREIGN KEY (`uniprot` )
    REFERENCES `uniprot` (`entry` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;
