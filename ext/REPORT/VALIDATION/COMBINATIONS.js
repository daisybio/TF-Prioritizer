let combinations = {COMBINATIONS}

function init() {
    let hm_selectors = document.getElementsByClassName("hm-selector");
    hm_selectors[0].click();

    select_first_possible_group();

    update_image();
}

function select_first_possible_group() {
    let group_selectors = document.getElementsByClassName("group-selector");
    let i;

    for (i = 0; i < group_selectors.length; i++) {
        if (!group_selectors[i].disabled) {
            group_selectors[i].click();
            break;
        }
    }
}

function hm_selector() {
    let hm_selectors = document.getElementsByClassName("hm-selector");
    let i;

    for (i = 0; i < hm_selectors.length; i++) {
        hm_selectors[i].addEventListener("click", function () {
            if (this.classList.contains("active")) return;

            let setClasses = !this.classList.contains("active");
            setClass(hm_selectors, "active", "remove");

            if (setClasses) {
                this.classList.toggle("active");
            }

            let possibleGroups = Object.keys(combinations[this.value]);

            let group_selectors = document.getElementsByClassName("group-selector");
            let j;

            for (j = 0; j < group_selectors.length; j++) {
                group_selectors[j].disabled = !possibleGroups.includes(group_selectors[j].value);
            }

            let activeGroup = get_active_element("group-selector");

            if (!(typeof activeGroup === 'undefined')) {
                if (!possibleGroups.includes(activeGroup.value)) {
                    select_first_possible_group();
                }
            }
        });
    }
}

function group_selector() {
    let group_selectors = document.getElementsByClassName("group-selector");
    let i;

    for (i = 0; i < group_selectors.length; i++) {
        group_selectors[i].addEventListener("click", function () {
            if (this.classList.contains("active")) return;

            let setClasses = !this.classList.contains("active");
            setClass(group_selectors, "active", "remove");

            if (setClasses) {
                this.classList.toggle("active");
            }


            let activeHm = get_active_element("hm-selector");

            let possibleGenes = combinations[activeHm.value][this.value];

            let dropdown = document.getElementById("select-gene");
            let l, k = dropdown.options.length - 1;

            for (l = k; l >= 0; l--) {
                dropdown.remove(l);
            }

            let m;
            for (m = 0; m < possibleGenes.length; m++) {
                let option = document.createElement("option");
                option.value = possibleGenes[m];
                option.textContent = possibleGenes[m];
                dropdown.appendChild(option);
            }

            update_image();
        });
    }
}

function update_image() {
    activeHm = get_active_element("hm-selector").value;
    activeGroup = get_active_element("group-selector").value;
    activeGene = document.getElementById("select-gene").value;

    let image = document.getElementById("validation-plot");
    image.src = activeHm + "/" + activeGroup + "/" + activeGene + ".png";
}

function get_active_element(className) {
    let elements = document.getElementsByClassName(className);
    let j;

    for (j = 0; j < elements.length; j++) {
        if (elements[j].classList.contains("active")) {
            return elements[j];
        }
    }
}