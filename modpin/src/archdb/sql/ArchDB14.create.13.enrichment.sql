DROP TABLE IF EXISTS `drugBank_class_enrichment` ;
DROP TABLE IF EXISTS `drugBank_subclass_enrichment` ;
DROP TABLE IF EXISTS `enzyme_subclass_enrichment` ;
DROP TABLE IF EXISTS `go_subclass_enrichment` ;
DROP TABLE IF EXISTS `scop_class_enrichment` ;
DROP TABLE IF EXISTS `enzyme_class_enrichment` ;
DROP TABLE IF EXISTS `scop_subclass_enrichment` ;
DROP TABLE IF EXISTS `go_class_enrichment` ;

CREATE  TABLE IF NOT EXISTS `drugBank_class_enrichment` (
  `class_nid` INT UNSIGNED NOT NULL ,
  `drugbank_id` VARCHAR(255) NOT NULL ,
  `enrichment_pvalue` VARCHAR(255) NULL ,
  `number_instances` INT NULL ,
  `frequency` VARCHAR(255) NULL ,
  `logodd` VARCHAR(255) NULL ,
  `mutual_information` VARCHAR(255) NULL ,
  INDEX `fk_drugBank_class_enrichment_cluster_class1_idx` (`class_nid` ASC) ,
  INDEX `fk_drugBank_class_enrichment_drugBank1_idx` (`drugbank_id` ASC) ,
  PRIMARY KEY (`class_nid`, `drugbank_id`) ,
  CONSTRAINT `fk_drugBank_class_enrichment_cluster_class1`
    FOREIGN KEY (`class_nid` )
    REFERENCES `cluster_class` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_drugBank_class_enrichment_drugBank1`
    FOREIGN KEY (`drugbank_id` )
    REFERENCES `drugBank` (`drugbank_id` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;

CREATE  TABLE IF NOT EXISTS `drugBank_subclass_enrichment` (
  `subclass_nid` INT UNSIGNED NOT NULL ,
  `drugbank_id` VARCHAR(255) NOT NULL ,
  `enrichment_pvalue` VARCHAR(255) NULL ,
  `number_instances` INT NULL ,
  `frequency` VARCHAR(255) NULL ,
  `logodd` VARCHAR(255) NULL ,
  `mutual_information` VARCHAR(255) NULL ,
  INDEX `fk_drugBank_subclass_enrichment_cluster_subclass1_idx` (`subclass_nid` ASC) ,
  INDEX `fk_drugBank_subclass_enrichment_drugBank1_idx` (`drugbank_id` ASC) ,
  PRIMARY KEY (`subclass_nid`, `drugbank_id`) ,
  CONSTRAINT `fk_drugBank_subclass_enrichment_cluster_subclass1`
    FOREIGN KEY (`subclass_nid` )
    REFERENCES `cluster_subclass` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_drugBank_subclass_enrichment_drugBank1`
    FOREIGN KEY (`drugbank_id` )
    REFERENCES `drugBank` (`drugbank_id` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;

CREATE  TABLE IF NOT EXISTS `enzyme_subclass_enrichment` (
  `subclass_nid` INT UNSIGNED NOT NULL ,
  `enzyme_id` VARCHAR(45) NOT NULL ,
  `enrichment_pvalue` VARCHAR(255) NULL ,
  `number_instances` INT NULL ,
  `frequency` VARCHAR(255) NULL ,
  `logodd` VARCHAR(255) NULL ,
  `mutual_information` VARCHAR(255) NULL ,
  INDEX `fk_drugBank_subclass_enrichment_cluster_subclass1_idx` (`subclass_nid` ASC) ,
  PRIMARY KEY (`subclass_nid`, `enzyme_id`) ,
  INDEX `fk_enzyme_subclass_enrichment_copy1_enzyme1_idx` (`enzyme_id` ASC) ,
  CONSTRAINT `fk_drugBank_subclass_enrichment_cluster_subclass10`
    FOREIGN KEY (`subclass_nid` )
    REFERENCES `cluster_subclass` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_enzyme_subclass_enrichment_copy1_enzyme1`
    FOREIGN KEY (`enzyme_id` )
    REFERENCES `enzyme` (`id` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;

CREATE  TABLE IF NOT EXISTS `go_subclass_enrichment` (
  `subclass_nid` INT UNSIGNED NOT NULL ,
  `GO_nid` INT UNSIGNED NOT NULL ,
  `enrichment_pvalue` VARCHAR(255) NULL ,
  `number_instances` INT NULL ,
  `frequency` VARCHAR(255) NULL ,
  `logodd` VARCHAR(255) NULL ,
  `mutual_information` VARCHAR(255) NULL ,
  INDEX `fk_drugBank_subclass_enrichment_cluster_subclass1_idx` (`subclass_nid` ASC) ,
  PRIMARY KEY (`subclass_nid`, `GO_nid`) ,
  INDEX `fk_GO_subclass_enrichment_copy1_GO1_idx` (`GO_nid` ASC) ,
  CONSTRAINT `fk_drugBank_subclass_enrichment_cluster_subclass11`
    FOREIGN KEY (`subclass_nid` )
    REFERENCES `cluster_subclass` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_GO_subclass_enrichment_copy1_GO1`
    FOREIGN KEY (`GO_nid` )
    REFERENCES `GO` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;

CREATE  TABLE IF NOT EXISTS `scop_subclass_enrichment` (
  `subclass_nid` INT UNSIGNED NOT NULL ,
  `scop_id` INT UNSIGNED NOT NULL ,
  `enrichment_pvalue` VARCHAR(255) NULL ,
  `number_instances` INT NULL ,
  `frequency` VARCHAR(255) NULL ,
  `logodd` VARCHAR(255) NULL ,
  `mutual_information` VARCHAR(255) NULL ,
  INDEX `fk_drugBank_subclass_enrichment_cluster_subclass1_idx` (`subclass_nid` ASC) ,
  PRIMARY KEY (`subclass_nid`, `scop_id`) ,
  INDEX `fk_SCOP_subclass_enrichment_scop_description1_idx` (`scop_id` ASC) ,
  CONSTRAINT `fk_drugBank_subclass_enrichment_cluster_subclass12`
    FOREIGN KEY (`subclass_nid` )
    REFERENCES `cluster_subclass` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_SCOP_subclass_enrichment_scop_description1`
    FOREIGN KEY (`scop_id` )
    REFERENCES `scop_description` (`id` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;

CREATE  TABLE IF NOT EXISTS `enzyme_class_enrichment` (
  `class_nid` INT UNSIGNED NOT NULL ,
  `enzyme_id` VARCHAR(45) NOT NULL ,
  `enrichment_pvalue` VARCHAR(255) NULL ,
  `number_instances` INT NULL ,
  `frequency` VARCHAR(255) NULL ,
  `logodd` VARCHAR(255) NULL ,
  `mutual_information` VARCHAR(255) NULL ,
  INDEX `fk_drugBank_class_enrichment_cluster_class1_idx` (`class_nid` ASC) ,
  PRIMARY KEY (`class_nid`, `enzyme_id`) ,
  INDEX `fk_enzyme_class_enrichment_enzyme1_idx` (`enzyme_id` ASC) ,
  CONSTRAINT `fk_drugBank_class_enrichment_cluster_class10`
    FOREIGN KEY (`class_nid` )
    REFERENCES `cluster_class` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_enzyme_class_enrichment_enzyme1`
    FOREIGN KEY (`enzyme_id` )
    REFERENCES `enzyme` (`id` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;

CREATE  TABLE IF NOT EXISTS `go_class_enrichment` (
  `class_nid` INT UNSIGNED NOT NULL ,
  `GO_nid` INT UNSIGNED NOT NULL ,
  `enrichment_pvalue` VARCHAR(255) NULL ,
  `number_instances` INT NULL ,
  `frequency` VARCHAR(255) NULL ,
  `logodd` VARCHAR(255) NULL ,
  `mutual_information` VARCHAR(255) NULL ,
  INDEX `fk_drugBank_class_enrichment_cluster_class1_idx` (`class_nid` ASC) ,
  PRIMARY KEY (`class_nid`, `GO_nid`) ,
  INDEX `fk_GO_class_enrichment_GO1_idx` (`GO_nid` ASC) ,
  CONSTRAINT `fk_drugBank_class_enrichment_cluster_class11`
    FOREIGN KEY (`class_nid` )
    REFERENCES `cluster_class` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_GO_class_enrichment_GO1`
    FOREIGN KEY (`GO_nid` )
    REFERENCES `GO` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;

CREATE  TABLE IF NOT EXISTS `scop_class_enrichment` (
  `class_nid` INT UNSIGNED NOT NULL ,
  `scop_id` INT UNSIGNED NOT NULL ,
  `enrichment_pvalue` VARCHAR(255) NULL ,
  `number_instances` INT NULL ,
  `frequency` VARCHAR(255) NULL ,
  `logodd` VARCHAR(255) NULL ,
  `mutual_information` VARCHAR(255) NULL ,
  INDEX `fk_drugBank_class_enrichment_cluster_class1_idx` (`class_nid` ASC) ,
  PRIMARY KEY (`class_nid`, `scop_id`) ,
  INDEX `fk_drugBank_class_enrichment_copy1_scop_description1_idx` (`scop_id` ASC) ,
  CONSTRAINT `fk_drugBank_class_enrichment_cluster_class12`
    FOREIGN KEY (`class_nid` )
    REFERENCES `cluster_class` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_drugBank_class_enrichment_copy1_scop_description1`
    FOREIGN KEY (`scop_id` )
    REFERENCES `scop_description` (`id` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;
