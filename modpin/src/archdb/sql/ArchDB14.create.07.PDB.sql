DROP function IF EXISTS `gethetero2PDB`;
DROP TABLE IF EXISTS `chain_contacts12` ;
DROP TABLE IF EXISTS `contacts12` ;
DROP TABLE IF EXISTS `chain_homology` ;
DROP TABLE IF EXISTS `hetero_contacts6` ;
DROP TABLE IF EXISTS `chain2uniprot` ;
DROP TABLE IF EXISTS `oldPDB` ;
DROP TABLE IF EXISTS `chain2taxid` ;
DROP TABLE IF EXISTS `chain2enzyme` ;
DROP TABLE IF EXISTS `chainrepeat` ;
DROP TABLE IF EXISTS `chainrepeatidx` ;
DROP TABLE IF EXISTS `site2position` ;
DROP TABLE IF EXISTS `site` ;
DROP TABLE IF EXISTS `hetero2PDB` ;
DROP TABLE IF EXISTS `chain` ;
DROP TABLE IF EXISTS `PDB` ;

CREATE  TABLE IF NOT EXISTS `PDB` (
  `pdb`        CHAR(4)      NOT NULL ,
  `date`       DATE         NOT NULL ,
  `header`     VARCHAR(255) NOT NULL ,
  `method`     ENUM('X-RAY DIFFRACTION','SOLUTION NMR','ELECTRON MICROSCOPY','NEUTRON DIFFRACTION','FIBER DIFFRACTION','ELECTRON CRYSTALLOGRAPHY','POWDER DIFFRACTION','SOLUTION SCATTERING','SOLID-STATE NMR','INFRARED SPECTROSCOPY','FLUORESCENCE TRANSFER') NOT NULL ,
  `resolution` FLOAT(3,2) NULL DEFAULT NULL COMMENT 'Value is NULL when NOT APPLICABLE' ,
  `Rfactor`    FLOAT(4,3) NULL DEFAULT NULL COMMENT 'Value is NULL when NOT APPLICABLE' ,
  `freeR`      FLOAT(4,3) NULL DEFAULT NULL COMMENT 'Value is NULL when NOT APPLICABLE' ,
  PRIMARY KEY (`pdb`) ,
  UNIQUE INDEX `pdb_id_UNIQUE` (`pdb` ASC) )
ENGINE = InnoDB;

CREATE  TABLE IF NOT EXISTS `chain` (
  `nid`   INT UNSIGNED   NOT NULL AUTO_INCREMENT ,
  `pdb`   CHAR(4)        NOT NULL ,
  `chain` CHAR(1) BINARY NOT NULL ,
  `name`  VARCHAR(255)       NULL DEFAULT NULL ,
  `type`  ENUM('P','N','O')  NULL DEFAULT 'O' ,
  `start` SMALLINT           NULL DEFAULT NULL ,
  `idxs`  CHAR(1)            NULL DEFAULT ' ' ,
  `end`   SMALLINT           NULL DEFAULT NULL ,
  `idxe`  CHAR(1)            NULL DEFAULT ' ' ,
  `clustered` ENUM('M','H','N') NULL DEFAULT 'N' ,
  INDEX `fk_chain_PDB_idx` (`pdb` ASC) ,
  PRIMARY KEY (`nid`) ,
  INDEX `pdb_chain_idx` (`pdb` ASC, `chain` ASC) ,
  INDEX `type_pdb_chain_idx` (`type` ASC, `pdb` ASC, `chain` ASC) ,
  UNIQUE INDEX `pdb_chain_uniqe` (`pdb` ASC, `chain` ASC) ,
  INDEX `type_nid` (`type` ASC, `nid` ASC) ,
  CONSTRAINT `fk_chain_PDB`
    FOREIGN KEY (`pdb` )
    REFERENCES `PDB` (`pdb` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;

CREATE  TABLE IF NOT EXISTS `hetero2PDB` (
  `nid` INT UNSIGNED NOT NULL AUTO_INCREMENT ,
  `chain` INT UNSIGNED NOT NULL ,
  `position` SMALLINT NOT NULL ,
  `hetero` VARCHAR(3) NOT NULL ,
  `inchain` TINYINT(1) NULL DEFAULT 0 ,
  INDEX `fk_hetero2PDB_hetero_idx` (`hetero` ASC) ,
  INDEX `fk_hetero2PDB_chain_idx` (`position` ASC) ,
  UNIQUE INDEX `id_UNIQUE` (`nid` ASC) ,
  PRIMARY KEY (`nid`) ,
  INDEX `fk_hetero2PDB_chain1_idx` (`chain` ASC) ,
  INDEX `inch_het` (`inchain` ASC, `hetero` ASC) ,
  INDEX `chain_inchain` (`chain` ASC, `inchain` ASC) ,
  UNIQUE INDEX `allpos_ind` (`chain` ASC, `position` ASC, `hetero` ASC) ,
  CONSTRAINT `fk_hetero2PDB_hetero`
    FOREIGN KEY (`hetero` )
    REFERENCES `hetero` (`id` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_hetero2PDB_chain1`
    FOREIGN KEY (`chain` )
    REFERENCES `chain` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;

CREATE  TABLE IF NOT EXISTS `site` (
  `nid` INT UNSIGNED NOT NULL AUTO_INCREMENT ,
  `pdb` CHAR(4) NOT NULL ,
  `name` VARCHAR(255) NOT NULL ,
  `description` VARCHAR(255) NULL DEFAULT NULL ,
  `bind` INT UNSIGNED NULL DEFAULT NULL ,
  PRIMARY KEY (`nid`) ,
  UNIQUE INDEX `nid_UNIQUE` (`nid` ASC) ,
  INDEX `fk_sites_PDB_idx` (`pdb` ASC) ,
  INDEX `fk_sites_hetero2PDB1_idx` (`bind` ASC) ,
  CONSTRAINT `fk_sites_PDB`
    FOREIGN KEY (`pdb` )
    REFERENCES `PDB` (`pdb` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_sites_hetero2PDB1`
    FOREIGN KEY (`bind` )
    REFERENCES `hetero2PDB` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;

CREATE  TABLE IF NOT EXISTS `site2position` (
  `site` INT UNSIGNED NOT NULL ,
  `chain` INT UNSIGNED NOT NULL ,
  `position` SMALLINT NOT NULL ,
  `idxp` CHAR(1) NOT NULL ,
  `type` VARCHAR(3) NOT NULL ,
  PRIMARY KEY (`site`, `position`, `chain`) ,
  INDEX `fk_site2position_chain2_idx` (`chain` ASC) ,
  INDEX `fk_site2position_hetero1_idx` (`type` ASC) ,
  CONSTRAINT `fk_site2position_site`
    FOREIGN KEY (`site` )
    REFERENCES `site` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_site2position_chain2`
    FOREIGN KEY (`chain` )
    REFERENCES `chain` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_site2position_hetero1`
    FOREIGN KEY (`type` )
    REFERENCES `hetero` (`id` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;

CREATE  TABLE IF NOT EXISTS `chainrepeatidx` (
  `nid` INT UNSIGNED NOT NULL AUTO_INCREMENT ,
  PRIMARY KEY (`nid`) )
ENGINE = InnoDB;

CREATE  TABLE IF NOT EXISTS `chainrepeat` (
  `groupid` INT UNSIGNED NOT NULL ,
  `chain` INT UNSIGNED NOT NULL ,
  INDEX `groupid_idx` (`groupid` ASC) ,
  INDEX `fk_chainrepeat_chain2_idx` (`chain` ASC) ,
  PRIMARY KEY (`groupid`, `chain`) ,
  CONSTRAINT `fk_chainrepeat_chainrepeatidx1`
    FOREIGN KEY (`groupid` )
    REFERENCES `chainrepeatidx` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_chainrepeat_chain2`
    FOREIGN KEY (`chain` )
    REFERENCES `chain` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;

CREATE  TABLE IF NOT EXISTS `chain2enzyme` (
  `chain` INT UNSIGNED NOT NULL ,
  `enzyme` VARCHAR(45) NOT NULL ,
  INDEX `fk_chain2enzyme_enzyme1_idx` (`enzyme` ASC) ,
  INDEX `fk_chain2enzyme_chain2_idx` (`chain` ASC) ,
  PRIMARY KEY (`chain`, `enzyme`) ,
  CONSTRAINT `fk_chain2enzyme_enzyme1`
    FOREIGN KEY (`enzyme` )
    REFERENCES `enzyme` (`id` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_chain2enzyme_chain2`
    FOREIGN KEY (`chain` )
    REFERENCES `chain` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;

CREATE  TABLE IF NOT EXISTS `chain2taxid` (
  `chain` INT UNSIGNED NOT NULL ,
  `taxid` INT UNSIGNED NOT NULL ,
  INDEX `fk_chain2taxid_taxid1_idx` (`taxid` ASC) ,
  INDEX `fk_chain2taxid_chain2_idx` (`chain` ASC) ,
  PRIMARY KEY (`chain`, `taxid`) ,
  CONSTRAINT `fk_chain2taxid_taxid1`
    FOREIGN KEY (`taxid` )
    REFERENCES `taxid` (`id` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_chain2taxid_chain2`
    FOREIGN KEY (`chain` )
    REFERENCES `chain` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;

CREATE  TABLE IF NOT EXISTS `oldPDB` (
  `oldid` CHAR(4) NOT NULL ,
  `pdb` CHAR(4) NOT NULL ,
  INDEX `fk_oldPDB_PDB1_idx` (`pdb` ASC) ,
  INDEX `odlPDB_old` (`oldid` ASC) ,
  CONSTRAINT `fk_oldPDB_PDB1`
    FOREIGN KEY (`pdb` )
    REFERENCES `PDB` (`pdb` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;

CREATE  TABLE IF NOT EXISTS `chain2uniprot` (
  `chain` INT UNSIGNED NOT NULL ,
  `uniprot` VARCHAR(255) NOT NULL ,
  `start` SMALLINT NULL DEFAULT NULL ,
  `idxs` CHAR(1) NULL DEFAULT ' ' ,
  `end` SMALLINT NULL DEFAULT NULL ,
  `idxe` CHAR(1) NULL DEFAULT ' ' ,
  INDEX `fk_chain2uniprot_chain2_idx` (`chain` ASC) ,
  INDEX `fk_chain2uniprot_uniprot1_idx` (`uniprot` ASC) ,
  INDEX `chain_unp_idx` (`chain` ASC, `uniprot` ASC) ,
  CONSTRAINT `fk_chain2uniprot_chain2`
    FOREIGN KEY (`chain` )
    REFERENCES `chain` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_chain2uniprot_uniprot1`
    FOREIGN KEY (`uniprot` )
    REFERENCES `uniprot` (`entry` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;

CREATE  TABLE IF NOT EXISTS `hetero_contacts6` (
  `nid` INT UNSIGNED NOT NULL AUTO_INCREMENT ,
  `chain` INT UNSIGNED NOT NULL ,
  `hetero` INT UNSIGNED NOT NULL ,
  `position` SMALLINT NOT NULL ,
  `idxp` CHAR(1) NOT NULL ,
  `type` VARCHAR(3) NOT NULL ,
  `distance` FLOAT(4,3) NOT NULL ,
  INDEX `fk_hetero_contacts6_hetero2PDB1_idx` (`hetero` ASC) ,
  PRIMARY KEY (`nid`) ,
  INDEX `fk_hetero_contacts6_chain2_idx` (`chain` ASC) ,
  INDEX `fk_hetero_contacts6_hetero1_idx` (`type` ASC) ,
  CONSTRAINT `fk_hetero_contacts6_hetero2PDB1`
    FOREIGN KEY (`hetero` )
    REFERENCES `hetero2PDB` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_hetero_contacts6_chain2`
    FOREIGN KEY (`chain` )
    REFERENCES `chain` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_hetero_contacts6_hetero1`
    FOREIGN KEY (`type` )
    REFERENCES `hetero` (`id` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;

CREATE  TABLE IF NOT EXISTS `chain_homology` (
  `master_nid` INT UNSIGNED NOT NULL ,
  `homolog_nid` INT UNSIGNED NOT NULL ,
  `master_length` INT UNSIGNED NULL ,
  `homolog_length` INT UNSIGNED NULL ,
  `homology` SMALLINT UNSIGNED NULL ,
  INDEX `fk_chain_homology_chain1_idx` (`master_nid` ASC) ,
  INDEX `fk_chain_homology_chain2_idx` (`homolog_nid` ASC) ,
  INDEX `master_homolog_idx` (`master_nid` ASC, `homolog_nid` ASC) ,
  CONSTRAINT `fk_chain_homology_chain1`
    FOREIGN KEY (`master_nid` )
    REFERENCES `chain` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_chain_homology_chain2`
    FOREIGN KEY (`homolog_nid` )
    REFERENCES `chain` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;

CREATE  TABLE IF NOT EXISTS `contacts12` (
  `nid` INT UNSIGNED NOT NULL AUTO_INCREMENT ,
  `alpha_distance` FLOAT(4,3) NOT NULL DEFAULT -1 ,
  `beta_distance` FLOAT(4,3) NOT NULL DEFAULT -1 ,
  `min_distance` FLOAT(4,3) NULL DEFAULT -1 ,
  `type` ENUM('PPI','PNI') NULL ,
  PRIMARY KEY (`nid`) ,
  INDEX `alpha` (`type` ASC, `alpha_distance` ASC, `nid` ASC) ,
  INDEX `beta` (`type` ASC, `beta_distance` ASC, `nid` ASC) ,
  INDEX `min` (`type` ASC, `min_distance` ASC, `nid` ASC) )
ENGINE = InnoDB;

CREATE  TABLE IF NOT EXISTS `chain_contacts12` (
  `contacts12` INT UNSIGNED NOT NULL ,
  `chain` INT UNSIGNED NOT NULL ,
  `position` SMALLINT NOT NULL ,
  `idxp` CHAR(1) NOT NULL ,
  `type` VARCHAR(3) NOT NULL ,
  INDEX `fk_chain_contacts12_contacts121_idx` (`contacts12` ASC) ,
  INDEX `fk_chain_contacts12_chain1_idx` (`chain` ASC) ,
  INDEX `fk_chain_contacts12_hetero1_idx` (`type` ASC) ,
  INDEX `contact_chain_res` (`contacts12` ASC, `chain` ASC, `position` ASC, `idxp` ASC) ,
  CONSTRAINT `fk_chain_contacts12_contacts121`
    FOREIGN KEY (`contacts12` )
    REFERENCES `contacts12` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_chain_contacts12_chain1`
    FOREIGN KEY (`chain` )
    REFERENCES `chain` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_chain_contacts12_hetero1`
    FOREIGN KEY (`type` )
    REFERENCES `hetero` (`id` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;

DELIMITER $$

CREATE FUNCTION`gethetero2PDB` (intchain INT, strposition SMALLINT, strhetero VARCHAR(255)) RETURNS INTEGER UNSIGNED
MODIFIES SQL DATA
BEGIN
  DECLARE identifier INTEGER UNSIGNED;
  DECLARE ciden INTEGER UNSIGNED;
  SELECT nid, COUNT(*) INTO identifier, ciden FROM hetero2PDB WHERE chain=intchain AND position=strposition AND hetero=strhetero;
  IF ciden = 0 THEN
    SELECT COUNT(*) INTO ciden FROM hetero WHERE id=strhetero;
    IF ciden > 0 THEN
      INSERT INTO hetero2PDB(chain,position,hetero) VALUES (intchain,strposition,strhetero);
      SET identifier = LAST_INSERT_ID();
    END IF;
  END IF;
  RETURN identifier;
END$$

DELIMITER ;

