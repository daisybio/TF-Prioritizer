import { ComponentFixture, TestBed } from '@angular/core/testing';

import { SingleTfAccordionComponent } from './single-tf-accordion.component';

describe('SingleTfAccordionComponent', () => {
  let component: SingleTfAccordionComponent;
  let fixture: ComponentFixture<SingleTfAccordionComponent>;

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      declarations: [ SingleTfAccordionComponent ]
    })
    .compileComponents();
  });

  beforeEach(() => {
    fixture = TestBed.createComponent(SingleTfAccordionComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
