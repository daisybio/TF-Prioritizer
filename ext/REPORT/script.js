function accordion() {
    let acc = document.getElementsByClassName("accordion");
    let i;

    for (i = 0; i < acc.length; i++) {
        acc[i].addEventListener("click", function () {
            this.classList.toggle("active");

            let panel = this.nextElementSibling;
            if (panel.style.display === "grid") {
                panel.style.display = "none";
            } else {
                panel.style.display = "grid";
            }
        });
    }
}

function selectImage(element, imageID) {
    let selectors = element.parentNode.children;

    let setClasses = !element.classList.contains("active");

    setClass(selectors, "active", "remove");

    let image = document.getElementById(imageID)

    image.src = element.value;

    if (setClasses) {
        element.classList.toggle("active");
    }
}

function openImageInTab(id) {
    let image = document.getElementById(id);
    window.open(image.src);
}

function openModal(modal_id) {
    let modal = document.getElementById(modal_id);
    modal.style.display = "block";
}

function closeModal(modal_id) {
    let modal = document.getElementById(modal_id);
    modal.style.display = "none";
}

function setClass(els, className, fnName) {
    for (let i = 0; i < els.length; i++) {
        els[i].classList[fnName](className);
    }
}

function tableMouseOver(gene, col, row) {
    colHead = document.getElementById(gene + "-col-" + col);
    rowHead = document.getElementById(gene + "-row-" + row);

    colHead.style.backgroundColor = "rgba(var(--black), 0.1)";
    rowHead.style.backgroundColor = "rgba(var(--black), 0.1)";
}

function tableMouseOut(gene, col, row) {
    colHead = document.getElementById(gene + "-col-" + col);
    rowHead = document.getElementById(gene + "-row-" + row);

    colHead.style.backgroundColor = "initial";
    rowHead.style.backgroundColor = "initial";
}

function init_selection(selection, combinations) {
    let group_selectors = document.getElementsByClassName(selection + " group-selector");
    select_group(selection, group_selectors[0], combinations);

    select_first_possible_subgroup(selection, combinations);

    update_image(selection, combinations);
}

function select_first_possible_subgroup(selection, combinations) {
    let subgroups = document.getElementsByClassName(selection + " subgroup-selector");
    let i;

    for (i = 0; i < subgroups.length; i++) {
        if (!subgroups[i].disabled) {
            select_subgroup(selection, subgroups[i], combinations);
            break;
        }
    }
}

function select_group(selection, element, combinations) {
    let groups = document.getElementsByClassName(selection + " group-selector");
    let subgroups = document.getElementsByClassName(selection + " subgroup-selector")

    if (element.classList.contains("active")) return;

    let inactive = !element.classList.contains("active");
    setClass(groups, "active", "remove");

    if (inactive) {
        element.classList.toggle("active");
    }

    let possible_subgroups;
    if (has_dropdown(selection)) {
        possible_subgroups = Object.keys(combinations[element.value]);
    } else {
        possible_subgroups = combinations[element.value];
    }


    let j;

    for (j = 0; j < subgroups.length; j++) {
        subgroups[j].disabled = !possible_subgroups.includes(subgroups[j].value);
    }

    let active_subgroup = get_active_element(selection + " subgroup-selector");

    if (!(typeof active_subgroup === 'undefined')) {
        if (!possible_subgroups.includes(active_subgroup.value)) {
            select_first_possible_subgroup(selection, combinations);
        } else {
            select_subgroup(selection, active_subgroup, combinations);
        }
    }
}

function select_subgroup(selection, element, combinations) {
    let group_selectors = document.getElementsByClassName(selection + " subgroup-selector");

    if (!element.classList.contains("active")) {
        let inactive = !element.classList.contains("active");
        setClass(group_selectors, "active", "remove");

        if (inactive) {
            element.classList.toggle("active");
        }
    }


    let active_group = get_active_element(selection + " group-selector");

    let possible_dropdown_values = combinations[active_group.value][element.value];

    if (has_dropdown(selection)) {
        let dropdown = document.getElementById(selection + "-dropdown");
        let l, k = dropdown.options.length - 1;

        for (l = k; l >= 0; l--) {
            dropdown.remove(l);
        }

        let m;
        for (m = 0; m < possible_dropdown_values.length; m++) {
            let option = document.createElement("option");
            option.value = possible_dropdown_values[m];

            let text;
            if (/^[0-9]+_*/.test(possible_dropdown_values[m])) {
                text = possible_dropdown_values[m].substring(possible_dropdown_values[m].split("_")[0].length + 1)
            } else {
                text = possible_dropdown_values[m];
            }

            option.textContent = (m + 1) + ". " + text;
            dropdown.appendChild(option);
        }

        let leftarrow = document.getElementById(selection + "-leftarrow");
        let rightarrow = document.getElementById(selection + "-rightarrow");

        let modalleftarrow = document.getElementById(selection + "-modal-leftarrow");
        let modalrightarrow = document.getElementById(selection + "-modal-rightarrow");

        dropdown.disabled = modalleftarrow.disabled = modalrightarrow.disabled = leftarrow.disabled = rightarrow.disabled = possible_dropdown_values.length < 2;
    }

    update_image(selection, combinations);
}

function move_dropdown(selection, delta, combinations) {
    if (has_dropdown(selection)) {
        let dropdown = document.getElementById(selection + "-dropdown");
        let options = dropdown.children;
        let i;

        for (i = 0; i < options.length; i++) {
            if (options[i].value === dropdown.value) {
                break;
            }
        }

        i += delta;

        if (i >= options.length) {
            i = 0;
        }
        if (i < 0) {
            i = options.length - 1;
        }

        dropdown.value = options[i].value;
    }

    update_image(selection, combinations);
}

function update_image(selection, combinations) {
    let active_group = get_active_element(selection + " group-selector").value;
    let active_subgroup = get_active_element(selection + " subgroup-selector").value;
    let file_name;
    let caption;

    if (has_dropdown(selection)) {
        let active_dropdown = document.getElementById(selection + "-dropdown").value;
        file_name = selection + "/" + active_group + "/" + active_subgroup + "/" + active_dropdown + ".png";
        let possibleGenes = combinations[active_group][active_subgroup];
        caption = (possibleGenes.indexOf(active_dropdown) + 1) + ". " + active_dropdown.split("_")[1];
    } else {
        file_name = selection + "/" + active_group + "/" + active_subgroup + ".png";
        caption = active_subgroup;
    }

    let image = document.getElementById(selection + "-image");
    let modal_image = document.getElementById(selection + "-modal-image");
    modal_image.src = image.src = file_name;

    let modal_caption = document.getElementById(selection + "-modal-caption");
    modal_caption.innerHTML = caption;
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

function open_all_panels() {
    let panels = document.getElementsByClassName("panel");

    let i;
    for (i = 0; i < panels.length; i++) {
        panels[i].style.display = "grid";
    }
}

function has_dropdown(selection) {
    return !!document.getElementById(selection + "-dropdown");
}