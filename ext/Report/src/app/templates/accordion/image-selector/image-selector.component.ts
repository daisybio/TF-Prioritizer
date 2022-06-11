import {Component, ElementRef, HostListener, Input, OnInit} from '@angular/core';
import {dataPrefix} from "../../../../assets/constants";

@Component({
  selector: 'image-selector',
  templateUrl: './image-selector.component.html',
  styleUrls: ['./image-selector.component.css']
})
export class ImageSelectorComponent implements OnInit {
  @Input()
  data!: {};

  @Input()
  visible: boolean = false;

  @Input()
  title: string = "undefined";

  @Input()
  sub: boolean = false;

  @Input()
  information: string | undefined;

  @Input()
  hideHeader: boolean = false;

  availableOptions: string[][] = [];
  activeOptions: string[] = [];
  allowedOptions: string[][] = [];

  hasDropDown: boolean = false;
  hasData: boolean = false;
  hasPlots: boolean = false;

  isEnabled: boolean = true;

  dropdownVisible: boolean = false;
  popupVisible: boolean = false;

  imageSource: string = "";

  constructor(private eRef: ElementRef) {
  }

  isAllowed(level: number, option: string) {
    return this.allowedOptions[level].includes(option);
  }

  ngOnInit(): void {
    this.collectLevels(this.data, 0);

    this.hasData = this.availableOptions[this.availableOptions.length - 1].includes("data");
    this.hasPlots = this.availableOptions[this.availableOptions.length - 1].includes("plot");

    this.isEnabled = this.hasData || this.hasPlots;

    if (this.isEnabled) {
      // remove last entries since they contain the "data" and "plot" entries
      this.availableOptions.pop();
      this.activeOptions.pop();

      if (this.activeOptions.length > 0) {
        this.hasDropDown = this.availableOptions[this.activeOptions.length - 1].length > 7;
        this.updateAllowedOptions(0);
      }
      this.updateImage();
    }
  }

  @HostListener('document:click', ['$event'])
  clickout(event: { target: any; }) {
    if (!this.eRef.nativeElement.contains(event.target)) {
      this.dropdownVisible = false;
    }
  }

  updateActiveOption(level: number, newValue: string) {
    this.activeOptions[level] = newValue;
    this.updateAllowedOptions(level + 1);
    this.updateImage();
  }

  updateImage() {
    let subMap = this.data;

    for (let level of this.activeOptions) {
      // @ts-ignore
      subMap = subMap[level];
    }

    // @ts-ignore
    this.imageSource = dataPrefix + subMap["plot"];
  }

  downloadData() {
    let subMap = this.data;

    for (let level of this.activeOptions) {
      // @ts-ignore
      subMap = subMap[level];
    }

    const link = document.createElement("a");
    link.setAttribute("target", "_blank");
    // @ts-ignore
    let file: string = subMap["data"];
    link.setAttribute("href", dataPrefix + file);
    link.setAttribute("download", "test.csv");
    document.body.appendChild(link);
    link.click();
    link.remove();
  }

  openNewTab(link: string) {
    window.open(link, '_blank');
  }

  updateAllowedOptions(level: number) {
    if (level >= this.activeOptions.length) {
      return;
    }

    this.allowedOptions[level] = [];

    for (let option of this.availableOptions[level]) {
      let subMap = this.getSubMap(level, option);
      if (this.containsDataLeafNodes(subMap)) {
        this.allowedOptions[level].push(option);
      }
    }

    if (!this.allowedOptions[level].includes(this.activeOptions[level])) {
      this.activeOptions[level] = this.allowedOptions[level][0];
    }

    if (level < this.activeOptions.length - 1) {
      this.updateAllowedOptions(level + 1);
    }
  }

  containsDataLeafNodes(map: {}) {
    if (!map) return false;
    if (this.isDataLeafNode(map)) return true;
    if (this.isLeafNode(map)) return false;

    for (let key of Object.keys(map)) {
      // @ts-ignore
      if (this.containsDataLeafNodes(map[key])) {
        return true;
      }
    }
    return false;
  }

  isDataLeafNode(map: {}) {
    return this.isLeafNode(map) && this.isDataNode(map);
  }

  isLeafNode(map: {}) {
    let keys = Object.keys(map);
    keys = keys.filter(entry => !["plot", "data"].includes(entry));

    return keys.length == 0;
  }

  isDataNode(map: {}) {
    return Object.keys(map).length > 0;
  }

  getSubMap(level: number, key: string) {
    let subMap = this.data;
    for (let i = 0; i < level; i++) {
      // @ts-ignore
      subMap = subMap[this.activeOptions[i]];
    }
    // @ts-ignore
    return subMap[key];
  }

  collectLevels(map: {}, level: number) {
    if (this.availableOptions.length <= level) {
      this.availableOptions[level] = [];
    }

    for (let key of Object.keys(map)) {
      if (!this.availableOptions[level].includes(key)) {
        this.availableOptions[level].push(key);
      }
      if (this.activeOptions.length <= level) {
        this.activeOptions.push(key);
      }

      if (!this.isLeafNode(map)) {
        // @ts-ignore
        let sub = map[key];
        this.collectLevels(sub, level + 1);
      }
    }
    this.availableOptions[level].sort();

    if (this.activeOptions.length <= level) {
      this.activeOptions.push(this.availableOptions[level][0]);
    }
  }

  moveDropdown(delta: number) {
    let dropDownLevel = this.activeOptions.length - 1;
    let currentIndex = this.allowedOptions[dropDownLevel].indexOf(this.activeOptions[dropDownLevel]);

    let newIndex = currentIndex + delta;

    if (newIndex < 0) {
      newIndex = this.allowedOptions[dropDownLevel].length - 1;
    } else if (newIndex >= this.allowedOptions[dropDownLevel].length) {
      newIndex = 0;
    }

    this.updateActiveOption(dropDownLevel, this.allowedOptions[dropDownLevel][newIndex]);
  }
}
