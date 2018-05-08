DROP TABLE IF EXISTS `PDBTM` ;
DROP TABLE IF EXISTS `PDBTM2chain` ;
DROP TABLE IF EXISTS `PDBTM_regions` ;

CREATE  TABLE IF NOT EXISTS `PDBTM_regions` (
  `id` CHAR NOT NULL ,
  `definition` VARCHAR(255) NULL ,
  PRIMARY KEY (`id`) )
ENGINE = InnoDB;

CREATE  TABLE IF NOT EXISTS `PDBTM` (
  `nid` INT NOT NULL AUTO_INCREMENT ,
  `tmres` FLOAT(6,3) NULL ,
  `type` VARCHAR(45) NULL ,
  `kwres` TINYINT(1) NULL ,
  PRIMARY KEY (`nid`) )
ENGINE = InnoDB;

CREATE  TABLE IF NOT EXISTS `PDBTM2chain` (
  `PDBTM` INT NOT NULL ,
  `chain` INT UNSIGNED NOT NULL ,
  `start` INT NOT NULL ,
  `idxs` CHAR(1) NULL ,
  `end` INT NULL ,
  `idxe` CHAR(1) NULL ,
  `region` CHAR NOT NULL ,
  INDEX `fk_PDBTM2chain_chain1_idx` (`chain` ASC) ,
  INDEX `fk_PDBTM2chain_PDBTM_regions1_idx` (`region` ASC) ,
  CONSTRAINT `fk_PDBTM2chain_PDBTM1`
    FOREIGN KEY (`PDBTM` )
    REFERENCES `PDBTM` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_PDBTM2chain_chain1`
    FOREIGN KEY (`chain` )
    REFERENCES `chain` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_PDBTM2chain_PDBTM_regions1`
    FOREIGN KEY (`region` )
    REFERENCES `PDBTM_regions` (`id` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;
