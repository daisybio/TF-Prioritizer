import {Component, Input, OnInit} from '@angular/core';

@Component({
  selector: 'image-selector',
  templateUrl: './image-selector.component.html',
  styleUrls: ['./image-selector.component.css']
})
export class ImageSelectorComponent implements OnInit {
  @Input()
    // @ts-ignore
  data: {};

  @Input()
  visible: boolean = false;

  availableOptions: string[][] = [];
  activeOptions: string[] = [];
  allowedOptions: string[][] = [];

  hasData: boolean = false;
  hasPlots: boolean = false;

  imageSource: string = "";
  dataPrefix: string = "assets/input/data"

  constructor() {
  }

  isAllowed(level: number, option: string) {
    return this.allowedOptions[level].includes(option);
  }

  ngOnInit(): void {
    this.collectLevels(this.data, 0);

    this.hasData = this.availableOptions[this.availableOptions.length - 1].includes("data");
    this.hasPlots = this.availableOptions[this.availableOptions.length - 1].includes("plot");

    this.availableOptions.pop();
    this.activeOptions.pop();

    this.updateAllowedOptions();
    this.updateImage();
  }

  updateActiveOption(level: number, newValue: string) {
    this.activeOptions[level] = newValue;
    this.updateAllowedOptions();
    this.updateImage();
  }

  updateAllowedOptions() {
    let subMap = this.data;
    let level = 0;

    for (let activeOption of this.activeOptions) {
      this.allowedOptions[level] = [];
      Object.keys(subMap).forEach(key => this.allowedOptions[level].push(key));

      if (!this.allowedOptions[level].includes(this.activeOptions[level])) {
        this.activeOptions[level] = this.allowedOptions[level][0];
      }

      level = level + 1;
      // @ts-ignore
      subMap = subMap[activeOption];
    }
  }

  updateImage() {
    let subMap = this.data;

    for (let level of this.activeOptions) {
      // @ts-ignore
      subMap = subMap[level];
    }

    // @ts-ignore
    this.imageSource = this.dataPrefix + subMap["plot"];
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
    link.setAttribute("href", this.dataPrefix + file);
    link.setAttribute("download", "test.csv");
    document.body.appendChild(link);
    link.click();
    link.remove();
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

      // @ts-ignore
      let sub = map[key];
      if (!(sub instanceof String || typeof sub == 'string')) {
        this.collectLevels(sub, level + 1);
      }
    }
  }
}
