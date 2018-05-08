DROP TABLE IF EXISTS `hetero_formula` ;
DROP TABLE IF EXISTS `hetero_parent` ;
DROP TABLE IF EXISTS `hetero` ;
DROP TABLE IF EXISTS `atom` ;

CREATE TABLE IF NOT EXISTS `atom` (
  `id`       VARCHAR(2)   NOT NULL ,
  `name`     VARCHAR(100)     NULL ,
  `atnumber` TINYINT UNSIGNED NULL ,
  PRIMARY KEY (`id`) )
ENGINE = InnoDB;

CREATE TABLE IF NOT EXISTS `hetero` (
  `id`      VARCHAR(3)   NOT NULL ,
  `name`    VARCHAR(255) NOT NULL ,
  `formula` VARCHAR(255) NOT NULL ,
  `weight`  FLOAT(6,3)   NOT NULL ,
  `charge`  SMALLINT     NOT NULL ,
  PRIMARY KEY (`id`) ,
  UNIQUE INDEX `id_UNIQUE` (`id` ASC) ,
  INDEX `charge_idx` (`charge` ASC, `id` ASC) )
ENGINE = InnoDB;

CREATE TABLE IF NOT EXISTS `hetero_parent` (
  `hetero` VARCHAR(3) NOT NULL ,
  `parent` VARCHAR(3) NOT NULL ,
  INDEX `fk_hetero_parent_hetero1_idx` (`hetero` ASC) ,
  INDEX `fk_hetero_parent_hetero2_idx` (`parent` ASC) ,
  PRIMARY KEY (`hetero`, `parent`) ,
  CONSTRAINT `fk_hetero_parent_hetero1`
    FOREIGN KEY (`hetero` )
    REFERENCES `hetero` (`id` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_hetero_parent_hetero2`
    FOREIGN KEY (`parent` )
    REFERENCES `hetero` (`id` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;

CREATE TABLE IF NOT EXISTS `hetero_formula` (
  `hetero` VARCHAR(3)   NOT NULL ,
  `atom`   VARCHAR(2)   NOT NULL ,
  `count`  TINYINT UNSIGNED NULL ,
  INDEX `fk_hetero_formula_hetero1_idx` (`hetero` ASC) ,
  INDEX `fk_hetero_formula_atom1_idx` (`atom` ASC) ,
  PRIMARY KEY (`hetero`, `atom`) ,
  CONSTRAINT `fk_hetero_formula_hetero1`
    FOREIGN KEY (`hetero` )
    REFERENCES `hetero` (`id` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_hetero_formula_atom1`
    FOREIGN KEY (`atom` )
    REFERENCES `atom` (`id` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;
