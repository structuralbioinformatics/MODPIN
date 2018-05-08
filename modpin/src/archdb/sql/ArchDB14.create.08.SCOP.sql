DROP TABLE IF EXISTS `chain2scop` ;
DROP TABLE IF EXISTS `scop` ;
DROP TABLE IF EXISTS `scop_description` ;

CREATE  TABLE IF NOT EXISTS `scop_description` (
  `id` INT UNSIGNED NOT NULL ,
  `type` ENUM('cl','cf','sf','fa','dm') NOT NULL ,
  `code` VARCHAR(255) NOT NULL ,
  `description` VARCHAR(255) NOT NULL ,
  PRIMARY KEY (`id`) ,
  UNIQUE INDEX `code_UNIQUE` (`id` ASC) )
ENGINE = InnoDB;

CREATE  TABLE IF NOT EXISTS `scop` (
  `domain` INT UNSIGNED NOT NULL ,
  `family` INT UNSIGNED NOT NULL ,
  `superfamily` INT UNSIGNED NOT NULL ,
  `fold` INT UNSIGNED NOT NULL ,
  `class` INT UNSIGNED NOT NULL ,
  PRIMARY KEY (`domain`) ,
  UNIQUE INDEX `domain_UNIQUE` (`domain` ASC) ,
  INDEX `fk_scop_scop_description2_idx` (`family` ASC) ,
  INDEX `fk_scop_scop_description3_idx` (`superfamily` ASC) ,
  INDEX `fk_scop_scop_description4_idx` (`fold` ASC) ,
  INDEX `fk_scop_scop_description5_idx` (`class` ASC) ,
  CONSTRAINT `fk_scop_scop_description1`
    FOREIGN KEY (`domain` )
    REFERENCES `scop_description` (`id` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_scop_scop_description2`
    FOREIGN KEY (`family` )
    REFERENCES `scop_description` (`id` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_scop_scop_description3`
    FOREIGN KEY (`superfamily` )
    REFERENCES `scop_description` (`id` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_scop_scop_description4`
    FOREIGN KEY (`fold` )
    REFERENCES `scop_description` (`id` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_scop_scop_description5`
    FOREIGN KEY (`class` )
    REFERENCES `scop_description` (`id` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;

CREATE  TABLE IF NOT EXISTS `chain2scop` (
  `chain` INT UNSIGNED NOT NULL ,
  `domain` INT UNSIGNED NOT NULL ,
  `start` SMALLINT NOT NULL,
  `idxs` CHAR(1) NULL DEFAULT ' ' ,
  `end` SMALLINT NULL DEFAULT NULL ,
  `idxe` CHAR(1) NULL DEFAULT ' ' ,
  `mapped` TINYINT(1) NOT NULL DEFAULT FALSE ,
  INDEX `fk_chain2scop_scop_idx` (`domain` ASC) ,
  INDEX `fk_chain2scop_chain1_idx` (`chain` ASC) ,
  CONSTRAINT `fk_chain2scop_scop`
    FOREIGN KEY (`domain` )
    REFERENCES `scop` (`domain` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_chain2scop_chain1`
    FOREIGN KEY (`chain` )
    REFERENCES `chain` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;


